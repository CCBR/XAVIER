import subprocess
from xavier.src.xavier.__main__ import main


def test_help():
    assert (
        "XAVIER"
        in subprocess.run(
            "./bin/xavier --help", capture_output=True, shell=True, text=True
        ).stdout
    )
