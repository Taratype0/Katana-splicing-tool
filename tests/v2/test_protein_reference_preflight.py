import hashlib
from types import SimpleNamespace

import pytest
from katana_next.domain.models import InputAsset
from katana_next.workflows.isoform_consequences import IsoformConsequencesWorkflow
from scripts.prepare_matched_protein_reference import reference_sequence_index


def test_reference_header_can_contain_arrow():
    raw=b'>sp|P1|test Enzyme A -> B\nMK\nK\n>sp|P2|test\nML\n'
    assert reference_sequence_index(raw)=={
        hashlib.sha256(b'MKK').hexdigest():['P1'], hashlib.sha256(b'ML').hexdigest():['P2']}
    with pytest.raises(ValueError,match='No UniProt'):
        reference_sequence_index(b'<html>Error</html>')


@pytest.mark.parametrize('problem', ['missing_species','wrong_species','wrong_hash'])
def test_protein_reference_rejected_in_plan_before_scan(tmp_path, problem):
    path=tmp_path/'protein.json'
    path.write_bytes(b'{}')
    asset=InputAsset(asset_id='protein',kind='uniprot_protein_json',path=path,
        metadata={'species':'Mus musculus','version':'1','source_url':'https://example.org/protein',
                  'sha256':'bad' if problem=='wrong_hash' else hashlib.sha256(path.read_bytes()).hexdigest()})
    species={} if problem=='missing_species' else {'species':'Homo sapiens' if problem=='wrong_species' else 'Mus musculus'}
    context=SimpleNamespace(manifest=SimpleNamespace(annotation=None,assets=[asset],metadata=species))
    with pytest.raises(ValueError,match='Declare project species|species mismatch|SHA256'):
        IsoformConsequencesWorkflow().input_paths(context,{})
