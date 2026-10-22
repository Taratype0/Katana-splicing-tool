"""Packaged-process API smoke test, with no native GUI automation claims."""
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen


def run_packaged_smoke(base_url: str, root: Path):
    import pandas as pd
    from katana_next.sashimi_fixture import create_sashimi_fixture
    if (root / "katana.project.json").exists():
        raise ValueError("Use a new empty fixture directory; never overwrite an acceptance project.")
    manifest = create_sashimi_fixture(root)

    def api(path, body=None):
        request = Request(base_url.rstrip("/")+path, data=json.dumps(body).encode() if body is not None else None,
                          headers={"Content-Type": "application/json"})
        with urlopen(request, timeout=60) as response:
            return json.load(response)

    project = api("/api/projects/open", {"manifest_path": str(manifest)})["snapshot"]["project_id"]
    prefix = f"/api/projects/{project}"
    jobs = []

    def run(workflow, params):
        job = api(prefix+f"/workflows/{workflow}/run", {"parameters": params})
        deadline = time.monotonic()+120
        while job["status"] in {"queued", "running"} and time.monotonic() < deadline:
            time.sleep(.2)
            job = api("/api/jobs/"+job["job_id"])
        if job["status"] not in {"complete", "cache_hit"}:
            raise RuntimeError(f"Packaged {workflow}: {job.get('error') or job['status']}")
        cached = api(prefix+f"/workflows/{workflow}/run", {"parameters": params})
        if cached["status"] != "cache_hit":
            raise RuntimeError("Packaged repeat did not hit cache")
        jobs.append({"workflow": workflow, "job": job, "repeat_cache_hit": True})
        return job

    run("single_cell.quantify", {"source": "bam", "min_cell_event_reads": 1})
    artifacts = api(prefix+"/artifacts")
    source = next(a for a in artifacts if a["name"] == "Single-cell junction choices")
    counts = pd.read_parquet(source["path"])
    if counts.inclusion.sum() != 18:
        raise RuntimeError("Packaged BAM molecule counts do not match the known fixture")
    selection = api(prefix+"/selections", {"name": "Packaged exact-junction smoke", "entity_type": "cell_splice_event",
                    "source_artifact_id": source["artifact_id"], "entity_ids": [counts.row_id.iloc[0]]})
    job = run("single_cell.sashimi", {"selection_id": selection["selection_id"], "flank": 20, "min_mapq": 20})
    artifacts = [a for a in api(prefix+"/artifacts") if a["artifact_id"] in job["artifact_ids"]]
    plots = [a for a in artifacts if a["metadata"].get("mime_type") == "image/png"]
    if len(plots) != 1:
        raise RuntimeError("Packaged Sashimi fixture should create exactly one donor page")
    with urlopen(base_url.rstrip("/")+prefix+f"/artifacts/{plots[0]['artifact_id']}/file") as response:
        if response.read(8) != b"\x89PNG\r\n\x1a\n":
            raise RuntimeError("Packaged preview endpoint did not return a PNG")
    all_jobs = api("/api/jobs?project_id="+project)
    if len(all_jobs) != 4:
        raise RuntimeError("Unexpected hidden jobs during packaged acceptance")
    result = {"project_id": project, "jobs": jobs, "png_endpoint_passed": True,
              "no_hidden_jobs": True, "scope": "Frozen-process API/BAM/PNG/cache smoke; synthetic data, not physical GUI acceptance"}
    (root / "packaged-smoke.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
