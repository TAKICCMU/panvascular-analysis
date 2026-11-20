"""This revision tests metadata only. Original full-preflight script is archived in audit/."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).with_name("test_metadata_revision.py")), run_name="__main__")
