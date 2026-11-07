"""Prepare a portable native Nextflow task, separate from the original run."""
from pathlib import Path
import json
import shutil
import zipfile
from .bundle import open_bundle


def prepare_task(request, output_dir, resources=None):
    output = Path(output_dir).resolve()
    if output.exists():
        raise ValueError('Use a new Nextflow task directory')
    with open_bundle(request) as bundle:
        request_id = bundle.request['request_id']
        requirements = bundle.request['resources']
        resource_ids = {r['resource_id'] for r in requirements}
        if set(resources or {}) - resource_ids:
            raise ValueError('Unknown resource binding')
        output.mkdir(parents=True)
        shutil.copytree(bundle.root, output / 'request')
    with zipfile.ZipFile(output / 'request.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted((output / 'request').rglob('*')):
            if path.is_file():
                archive.write(path, path.relative_to(output / 'request').as_posix())
    # This is the very same module used by pipeline/main.nf -entry TARGETS.
    checkout = Path(__file__).resolve().parents[2]
    module = checkout / 'pipeline/workflows/katana_targets.nf'
    if not module.is_file():
        module = Path(__file__).resolve().parents[1] / 'flow/_pipeline/workflows/katana_targets.nf'
    shutil.copy2(module, output / 'katana_targets.nf')
    main = (checkout / 'pipeline/main.nf')
    if not main.is_file():
        main = Path(__file__).resolve().parents[1] / 'flow/_pipeline/main.nf'
    entry = main.read_text(encoding='utf-8').split('// A new follow-up run from a portable request exported by the browser.')[1]
    entry = entry.replace("'./workflows/katana_targets'", "'./katana_targets'")
    (output / 'main.nf').write_text(entry, encoding='utf-8')
    (output / 'nextflow.config').write_text("nextflow.enable.dsl = 2\nparams.katana_container = null\nparams.target_request = 'request.zip'\nparams.target_bindings = null\nparams.outdir = 'results'\ntrace.enabled = true\ntrace.file = 'trace.txt'\n", encoding='utf-8')
    template = {'schema':'katana.target_bindings/v1', 'resources':{r['resource_id']: '' for r in requirements}}
    (output / 'bindings.template.json').write_text(json.dumps(template, indent=2), encoding='utf-8')
    (output / 'resource-requirements.json').write_text(json.dumps(requirements, indent=2), encoding='utf-8')
    # Execution-machine paths stay outside the portable task archive.
    if resources:
        (output / 'bindings.local.json').write_text(json.dumps({'schema':template['schema'],'resources':resources}, indent=2), encoding='utf-8')
    command = 'nextflow run main.nf -entry TARGETS --target_request request.zip --target_bindings bindings.json --outdir results -work-dir work'
    (output / 'README.txt').write_text('Katana selected-target follow-up task\n\nInstall this Katana version on the execution host. Fill bindings.template.json with execution-host BAM/index/GTF paths and save as bindings.json. All files are staged as native Nextflow path inputs.\n\n' + command + '\n\nUse -c site.config for your executor/resources and -resume to resume. No upstream analysis is launched. Each target is published under results/targets/<target_id>, with independent completion/partial/failure records, REPORT.zip and result-manifest.json. Reopen each returned result directory in the Gene / Event workspace. Missing read resources remain unavailable.\n', encoding='utf-8')
    archive_path = output / 'NEXTFLOW_TASK.zip'
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in ('request.zip','main.nf','katana_targets.nf','nextflow.config','bindings.template.json','resource-requirements.json','README.txt'):
            archive.write(output / name, name)
    result = {'status':'prepared', 'execution':'native_nextflow', 'request_id':request_id,
              'directory':str(output), 'archive':str(archive_path), 'command':command,
              'submitted':False}
    (output / 'task.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result
