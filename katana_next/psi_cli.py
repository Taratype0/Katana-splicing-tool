"""Fit PSI networks and project frozen modules using the same production code as the app."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import pandas as pd

from katana_cosplice.config import NetworkConfig
from katana_cosplice.network import fit_network
from katana_next.science.psi_network import PsiFilterConfig, psi_filter_audit, prepare_psi_network, fit_psi_network,fit_psi_backbone
from katana_next.science.frozen_psi import freeze_psi_modules, save_frozen_model, load_frozen_model, project_psi_modules, _samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    for command in ("fit", "project"):
        p = subs.add_parser(command)
        for flag in ("psi", "samples", "out", "feature-space-id"):
            p.add_argument("--"+flag, required=True)
        p.add_argument("--traits", help="Numeric sample-by-trait TSV, first column sample_id")
        p.add_argument('--reference-scores',help='Optional numeric sample-by-reference TSV; use the same frozen reference transform for projection')
        p.add_argument('--reference-space-id',help='Identity of the frozen reference transform/feature definitions')
        if command == "fit":
            p.add_argument("--trait-prediction-method", choices=("single_module", "multimodule_ridge",'validated_hierarchy'), default="validated_hierarchy")
            p.add_argument('--filter-profile',choices=('legacy','coverage_variable'),default='coverage_variable')
            p.add_argument('--max-features',type=int,default=2000)
            p.add_argument("--rscript")
            p.add_argument("--method", choices=("bicor", "pearson"), default="bicor")
            p.add_argument("--unsigned", action="store_true")
            p.add_argument("--soft-power", type=int, default=6)
            p.add_argument("--backbone-only", action="store_true", help="Explicitly skip the 50 confidence refits; does not alter membership")
            for name, value in asdict(PsiFilterConfig()).items():
                p.add_argument("--"+name.replace("_", "-"), type=float, default=value)
        else:
            p.add_argument("--projection-method", choices=('auto',"mean_impute", "observed_loading_ls",'guarded_loading'), default='auto')
            p.add_argument("--min-loading-fraction", type=float, default=.8)
            p.add_argument("--model", required=True)
            p.add_argument("--min-observed-fraction", type=float, default=.8)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    matrix = pd.read_csv(args.psi, sep="\t", index_col=0)
    samples = pd.read_csv(args.samples, sep="\t")
    traits = pd.read_csv(args.traits, sep="\t", index_col=0) if args.traits else None
    reference=pd.read_csv(args.reference_scores,sep='\t',index_col=0) if args.reference_scores else None
    _samples(matrix, samples)
    if args.command == "fit":
        if reference is not None and (traits is None or args.trait_prediction_method!='validated_hierarchy'):
            raise ValueError('Reference fallback requires training traits and validated_hierarchy')
        filters = PsiFilterConfig(**{k: getattr(args, k) for k in asdict(PsiFilterConfig())})
        qc = psi_filter_audit(matrix, filters)
        if args.filter_profile=='coverage_variable':
            from katana_next.science.psi_robust import coverage_variable_audit
            qc=coverage_variable_audit(matrix,filters,max_features=args.max_features)
        qc.to_csv(out/"event_filter_audit.tsv", sep="\t", index_label="event_key")
        network = NetworkConfig(wgcna_rscript=args.rscript, method=args.method, signed=not args.unsigned, soft_power=args.soft_power)
        if not 1 <= args.soft_power <= 30:
            raise ValueError("soft_power must be within [1,30]")
        if args.backbone_only:
            measured, complete, qc, means = prepare_psi_network(matrix, samples, filters=filters,
                filter_profile=args.filter_profile,max_features=args.max_features)
            baseline = fit_psi_backbone(complete,qc,network,samples)
        else:
            measured, qc, means, baseline, confidence = fit_psi_network(matrix, samples, network=network, filters=filters,
                progress=print,filter_profile=args.filter_profile,max_features=args.max_features)
            confidence.assignments.to_csv(out/"confidence_assignments.tsv", sep="\t", index=False)
            confidence.bootstrap_runs.to_csv(out/"confidence_refits.tsv", sep="\t", index=False)
            confidence.module_summary.to_csv(out/"confidence_module_summary.tsv", sep="\t", index=False)
        model = freeze_psi_modules(measured, means, baseline, samples, feature_space_id=args.feature_space_id,
                                  network=network, filters=filters, traits=traits,
                                  trait_prediction_method=args.trait_prediction_method if traits is not None else 'single_module',
                                  reference_scores=reference,reference_space_id=args.reference_space_id,filter_profile=args.filter_profile)
        save_frozen_model(model, out/"frozen_psi_modules.json")
        baseline.assignments.to_csv(out/"module_assignments.tsv", sep="\t", index=False)
        baseline.program_scores.to_csv(out/"discovery_scores.tsv", sep="\t")
        baseline.edges.to_csv(out/"network_edges.tsv", sep="\t", index=False)
        baseline.diagnostics.to_csv(out/"network_diagnostics.tsv", sep="\t", index=False)
        measured.to_csv(out/"measured_discovery_psi.tsv.gz", sep="\t")
        (out/"resolved_configuration.json").write_text(json.dumps({"network": asdict(network), "filters": asdict(filters),
            'filter_profile':args.filter_profile,'max_features':args.max_features,'reference_space_id':args.reference_space_id}, indent=2), encoding="utf-8")
        summary = dict(status="complete", input_events=len(matrix), retained_events=len(measured),
            modules=len(model["modules"]),network_status=model['network_status'],
            confidence_refits=0 if args.backbone_only or not model['modules'] else 50, model_sha256=model["sha256"])
    else:
        model = load_frozen_model(args.model)
        projection=args.projection_method
        if projection=='auto':projection='guarded_loading' if model.get('schema')=='katana-frozen-psi-1.2' else 'mean_impute'
        outputs = project_psi_modules(model, matrix, samples, feature_space_id=args.feature_space_id,
                                     min_observed_fraction=args.min_observed_fraction, traits=traits,
                                     projection_method=projection, min_loading_fraction=args.min_loading_fraction,
                                     reference_scores=reference,reference_space_id=args.reference_space_id)
        for name, frame in outputs.items():
            frame.to_csv(out/(name+".tsv"), sep="\t", index=name == "scores")
        save_frozen_model(model, out/"frozen_psi_modules.json")
        summary = dict(status="complete", validation_samples=len(matrix.columns), modules=len(model["modules"]),
                       module_refit=False, model_sha256=model["sha256"], regulatory_accuracy="not_tested")
    (out/"COMPLETE.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
