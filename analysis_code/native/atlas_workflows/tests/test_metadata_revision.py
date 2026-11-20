"""Run with Python standard library. No h5ad or server expressions are changed."""
from pathlib import Path
import csv, json, hashlib
B = Path(__file__).resolve().parents[1]
with (B / "inputs/metadata_h5ad_patch_v4.csv").open(encoding="utf-8-sig", newline="") as f:
    m = list(csv.DictReader(f))
assert len(m) == 247
assert len({r["sample_uid"] for r in m}) == 247
assert all(r["donor_uid"] for r in m)
assert len({r["donor_uid"] for r in m}) == 168
assert sum(not r["age"] for r in m) == 76
assert sum(not r["sex"] for r in m) == 56
for study, count in [("GSE166676", 6), ("GSE184073", 2)]:
    s = [r for r in m if r["study"] == study]
    assert len(s) == count and len({r["donor_uid"] for r in s}) == count
    assert all(r["donor_uid"] == study + "::DONOR::" + r["source_sample_accession"] for r in s)
q = [r for r in m if r["study"] == "GSE237230"]
assert len(q) == 8 and len({r["donor_uid"] for r in q}) == 4
for p in [("GSM7597522","GSM7597523"), ("GSM7597524","GSM7597525"),
          ("GSM7597526","GSM7597527"), ("GSM7597528","GSM7597529")]:
    assert len({r["donor_uid"] for r in q if r["source_sample_accession"] in p}) == 1
for rel, expected in json.loads((B / "BUNDLE_CHECKSUMS.json").read_text()).items():
    actual = hashlib.sha256((B / rel).read_bytes()).hexdigest()
    assert actual == expected, rel
print("PASS: 247 samples; no missing donor keys; 168 keys; only confirmed 8-sample donor update; checksum manifest valid.")
