import subprocess
import sys
from pathlib import Path


def test_run_script_help_works() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(repo_root / "pipeline" / "run.py"), "--help"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "pipeline stage" in result.stdout.lower()
