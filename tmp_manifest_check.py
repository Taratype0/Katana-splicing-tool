import traceback
try:
    from src.services.project_service import ProjectService
    print('imported service')
    svc = ProjectService()
    print('created service')
    svc.load_project(r'E:\Project\RNASEQ\Result\Ythdc1_backup', force_auto_scan=True)
    print('loaded project')
    svc.confirm_project()
    svc.confirm_pairing()
    svc.confirm_comparison_sets()
    print('confirmed flow')
    svc.run_main_analysis_with_progress()
    print('ran analysis')
    manifest = svc.get_or_build_sashimi_manifest()
    print('manifest_rows', len(manifest))
    if not manifest.empty:
        print(manifest.head().to_string())
    else:
        print('EMPTY')
except Exception:
    traceback.print_exc()
    raise
