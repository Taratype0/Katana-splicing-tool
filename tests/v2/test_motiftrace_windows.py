import base64
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_windows import FIELDS, METHODS, Q, MatchedWindowBundle, render_windows, support_status
from katana_next.science.rbp_sequence_results import digest


def fixture(tmp_path):
    root = tmp_path / "saved"
    root.mkdir()
    events, cases, members = [], {"PBS": {}}, {}
    for i, role in enumerate(("up", "dn", "bg")):
        coords = dict(exonStart_0base=100+i, exonEnd=130+i, upstreamES=50, upstreamEE=80, downstreamES=170, downstreamEE=200)
        token = "SE|chr=chr1|strand=+|"+"|".join(f"{k}={v}" for k,v in coords.items())
        pd.DataFrame([dict(event_id=token, comparison_id="PBS", gene_id=role, chr="chr1", strand="+", dpsi=.2, fdr=.01, **coords)]).to_parquet(root/(role+".parquet"))
        cases["PBS"][role] = role
        members[role] = role+".parquet"
        events.append(dict(event_id=token, gene_id=role, chr="chr1", strand="+", role=role, dpsi=.2, DAS_FDR=.01))
    rows = []
    for direction in ("inclusion_up", "inclusion_down"):
        for method in METHODS:
            for position in (0,1):
                rows.append(dict(contrast="PBS", motif_index=98, RBP="PTBP1", native_RBP="PTBP1", pattern="CTCT[CT][CT]",
                    region_index=2, region="UpstreamIntron", native_position=position,
                    window_start_RNA_relative=-100+position, window_end_RNA_relative=-50+position,
                    direction=direction, arm=method, raw_P=.001, **{Q:.01 if method == METHODS[0] else .2},
                    primary_family_P=.001, n_foreground=30, n_background=31, mean_foreground=.2, mean_background=.1, mean_difference=.1))
    pd.DataFrame(rows).to_parquet(root/"tests.parquet")
    plan = dict(threshold=.05, author_source_modified=False, revision="fixture_not_real",
        parameters=dict(intron=200,exon=50,window=50,step=1,stat_method="mannwhitney_greater"),
        analysis_arms={}, FDR_primary="separate", canonical_masks="unchanged", validity="NA retained", limitations=[])
    (root/"plan.json").write_text(json.dumps(plan))
    pd.DataFrame([dict(window_id="2:0:inclusion_down:98", dataset="primary", quality="author", RBP="PTBP1", pattern="CTCT[CT][CT]",
        relative_start=-100,relative_stop=-50,original_author_P=.001,original_author_q=.01,original_MotifTrace_P=.001,original_MotifTrace_q=.2,
        foreground_motif_CLIP=2, foreground_CLIP=3,background_CLIP=1, foreground_N=30,background_N=31,
        external_CLIP_raw_P=.001,external_CLIP_q_all_available_windows=.02)]).to_csv(root/"external.tsv",sep="\t",index=False)
    members.update(tests="tests.parquet",plan="plan.json",external_windows="external.tsv")
    reg = dict(schema="motiftrace.matched_windows.v1", state="complete",species="Mus_musculus",genome_build="GRCm39",
        source_bundle_manifest_sha256="source",threshold=.05,families=dict.fromkeys(METHODS,4),cases=cases,
        external_contrast="PBS",external_sources={"primary":{"level":"primary_author_significant","quality":"author"}},members={})
    def seal():
        reg["members"] = {k:dict(path=p,bytes=(root/p).stat().st_size,sha256=digest(root/p)) for k,p in members.items()}
        (root/"registration.json").write_text(json.dumps(reg))
        return SimpleNamespace(path=root/"registration.json", metadata=dict(bundle_root=str(root),sha256=digest(root/"registration.json")))
    return root, reg, {"cases":["PBS"],"events":{"PBS":events}}, seal


def load(asset):
    return MatchedWindowBundle(asset,species="Mus_musculus",genome_build="GRCm39",source_manifest_sha256="source")


def test_original_method_values_are_binary_exact_and_support_is_not_significance(tmp_path):
    _,_,base,seal=fixture(tmp_path)
    result=load(seal()).payload(base)
    values=np.frombuffer(base64.b64decode(result["data"]["PBS|98|inclusion_down|Katana_overlapping_coverage"]),dtype="<f8").reshape(-1,len(FIELDS))
    assert values[0,:2].tolist()==[.001,.2]
    assert result["families"]==dict.fromkeys(METHODS,4)
    assert result["evidence"]["PBS|2:0:inclusion_down:98"]["status"]==["reference_overlap_method_passed","reference_enriched_method_not_passed"]
    assert result["new_inference"] is False


@pytest.mark.parametrize("mutation",["species","source","hash","family","coordinate","effect","background","probability","original_q","grid","duplicate"])
def test_reject_changed_binding_or_probabilities(tmp_path,mutation):
    root,reg,base,seal=fixture(tmp_path)
    if mutation=="species":reg["species"]="Homo_sapiens"
    elif mutation=="source":reg["source_bundle_manifest_sha256"]="other"
    elif mutation=="family":reg["families"][METHODS[0]]=1
    elif mutation=="coordinate":
        d=pd.read_parquet(root/"up.parquet");d.loc[0,"exonEnd"]+=1;d.to_parquet(root/"up.parquet")
    elif mutation=="effect":base["events"]["PBS"][0]["dpsi"]=-.2
    elif mutation=="background":base["events"]["PBS"].pop()
    elif mutation in ("probability","grid","duplicate"):
        d=pd.read_parquet(root/"tests.parquet")
        if mutation=="probability":d.loc[0,Q]=2
        elif mutation=="grid":d=d.iloc[1:]
        else:d=pd.concat([d,d.iloc[:1]])
        d.to_parquet(root/"tests.parquet")
    elif mutation=="original_q":
        d=pd.read_csv(root/"external.tsv",sep="\t");d.loc[0,"original_MotifTrace_q"]=.04;d.to_csv(root/"external.tsv",sep="\t",index=False)
    asset=seal()
    if mutation=="hash":(root/"tests.parquet").write_bytes(b"tampered")
    with pytest.raises((ValueError, KeyError)):
        load(asset).payload(base)


@pytest.mark.parametrize("level,q,expected",[("ambiguous_family",.9,"family_reference_only"),("primary_author_significant",None,"reference_overlap_method_NA"),("secondary_reproducible",.9,"reference_overlap_method_not_passed")])
def test_evidence_levels_NA_and_secondary_are_not_upgraded(level,q,expected):
    evidence=dict(dataset="x",level=level,motif_overlap=True,external_q=.001)
    assert support_status(evidence,q)==expected
    assert support_status(None,.9)=="not_assessed"
    assert support_status(dict(dataset=None),.9)=="no_qualified_reference"


def test_portable_report_escapes_injected_metadata_and_has_no_network_worker(tmp_path):
    data={"warning":"</script><img src=x onerror=BAD>"}
    page=render_windows(data,tmp_path/"viewer.html").read_text(encoding="utf-8")
    assert "</script><img" not in page
    assert "\\u003c/script>" in page
    assert "fetch(" not in page and "XMLHttpRequest" not in page
