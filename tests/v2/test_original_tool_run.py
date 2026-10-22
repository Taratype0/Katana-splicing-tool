from pathlib import Path
import subprocess
import sys

import pytest

from katana_next.infrastructure.original_tool_run import run_original_tool


def repository(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", str(repo)], check=True, capture_output=True)
    (repo / "source.txt").write_text("original")
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.org", "commit", "-m", "fixture"], cwd=repo, check=True, capture_output=True)
    return repo


def test_verified_outputs_and_immutable_run(tmp_path):
    repo = repository(tmp_path)
    out = tmp_path / "result"
    args = dict(output=out, repository=repo, inputs=[repo / "source.txt"],
                commands=[[sys.executable, "-c", f"from pathlib import Path; Path({str(out / 'answer.tsv')!r}).write_text('x\\n1\\n')"]],
                expected_outputs=["answer.tsv"], boundary="test")
    result = run_original_tool(**args)
    assert result["status"] == "execution_complete"
    assert result["scientific_validation"] == "not_asserted_by_executor"
    assert (out / "completed.json").is_file()
    with pytest.raises(ValueError, match="fresh"):
        run_original_tool(**args)


def test_missing_output_is_failure(tmp_path):
    repo = repository(tmp_path)
    out = tmp_path / "result"
    with pytest.raises(FileNotFoundError):
        run_original_tool(output=out, repository=repo, inputs=[],
            commands=[[sys.executable, "-c", "pass"]], expected_outputs=["missing"], boundary="test")
    assert (out / "failed.json").exists()
    assert not (out / "completed.json").exists()


def test_original_source_and_output_boundary(tmp_path):
    repo = repository(tmp_path)
    args = dict(output=tmp_path / "result", repository=repo, inputs=[],
                commands=[[sys.executable, "-c", "pass"]], expected_outputs=["../escape"], boundary="test")
    with pytest.raises(ValueError, match="inside"):
        run_original_tool(**args)
    args["expected_outputs"] = ["result.tsv"]
    (repo / "source.txt").write_text("patched")
    with pytest.raises(ValueError, match="modifications"):
        run_original_tool(**args)
