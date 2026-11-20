from pathlib import Path
import importlib.util

L2_UPDATE_ROOT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908"
)
L2_UPDATE_SCRIPT = L2_UPDATE_ROOT / "L2_followup_metadata_v4_code/run_L2_followup_and_update.py"
if not L2_UPDATE_SCRIPT.is_file():
    raise FileNotFoundError(L2_UPDATE_SCRIPT)
_v4_spec = importlib.util.spec_from_file_location("l2_v4_donor_revision", L2_UPDATE_SCRIPT)
l2_v4_donor_revision = importlib.util.module_from_spec(_v4_spec)
_v4_spec.loader.exec_module(l2_v4_donor_revision)

# Prefer the last successfully returned complete h5ad if present.
# The original driver rejects conflicts with newer, incompatible L2 labels.
V4_DONOR_SOURCE = globals().get(
    "L2_UPDATED_H5AD",
    L2_UPDATE_ROOT / "panvascular_CPU_theta2_level1_frozen_v1.h5ad",
)
V4_DONOR_RUN, V4_DONOR_H5AD = l2_v4_donor_revision.run_from_notebook(
    root=L2_UPDATE_ROOT,
    source=Path(V4_DONOR_SOURCE),
    macrophage_run=globals().get("MACROPHAGE_L2_RUN"),
    other_run=globals().get("OTHER_L2_RUN"),
    total_threads=40,
    max_workers=2,
    resume_run=None,       # Do not reuse an old metadata snapshot.
    run_targeted=False,   # This invocation updates full-h5ad metadata only.
)
print("Updated full h5ad:", V4_DONOR_H5AD)
