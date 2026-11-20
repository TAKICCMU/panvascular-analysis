#!/usr/bin/env python3
"""Recover timed-out GTEx requests from the official imported GTEx v8 files.

The uniformly processed eQTL Catalogue files are primary.  This script is a
source-documented fallback only for requests marked TABIX_FAILED in the primary
manifest.  KCNK3 files that were already downloaded from the same imported
release are reused after content validation; other failed contexts are queried
from the corresponding complete tissue file with remote tabix.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from pathlib import Path


ROOT = Path("/path/to/panvas_work/Tier1_full_colocalisation_v1")
MANIFEST = ROOT / "00_metadata" / "qtl_download_manifest.tsv"
OUT = ROOT / "01_qtl_regions_uniform"
OLD_KCNK3 = Path("/path/to/panvas_work/GWAS_eQTL_colocalisation/01_qtl_regions")
TABIX = "/opt/homebrew/bin/tabix"
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge"

COLS = [
    "variant", "r2", "pvalue", "molecular_trait_object_id",
    "molecular_trait_id", "maf", "gene_id", "median_tpm", "beta", "se",
    "an", "ac", "chromosome", "position", "ref", "alt", "type", "rsid",
]

CONTEXT_FILE = {
    "Artery_Aorta": "Artery_Aorta.tsv.gz",
    "Artery_Coronary": "Artery_Coronary.tsv.gz",
    "Artery_Tibial": "Artery_Tibial.tsv.gz",
    "Cells_Cultured_fibroblasts": "Cells_Cultured_fibroblasts.tsv.gz",
    "Heart_Atrial_Appendage": "Heart_Atrial_Appendage.tsv.gz",
    "Heart_Left_Ventricle": "Heart_Left_Ventricle.tsv.gz",
    "Whole_Blood": "Whole_Blood.tsv.gz",
}


def validate(path: Path, gene_id: str) -> int:
    if not path.exists() or path.stat().st_size <= 70:
        return 0
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != COLS:
            raise ValueError(f"Unexpected columns in {path}: {reader.fieldnames}")
        return sum(row["gene_id"].split(".", 1)[0] == gene_id for row in reader)


def fetch(row: dict[str, str]) -> tuple[Path, int, str]:
    gene = row["gene"]
    context = row["context"]
    gene_id = row["gene_id"]
    target = OUT / f"{gene}__{context}__cis1Mb.tsv.gz"
    done = Path(str(target) + ".done")

    old = OLD_KCNK3 / target.name
    old_done = Path(str(old) + ".done")
    if gene == "KCNK3" and old.exists() and old_done.exists():
        n = validate(old, gene_id)
        if n:
            shutil.copy2(old, target)
            done.write_text(f"complete\t{n}\n")
            return target, n, f"{BASE}/{CONTEXT_FILE[context]}"

    chrom = row["region_query"].split(":", 1)[0]
    region = row["region_query"]
    url = f"{BASE}/{CONTEXT_FILE[context]}"
    proc = subprocess.Popen(
        [TABIX, url, region], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, bufsize=1,
    )
    assert proc.stdout is not None
    tmp = Path(str(target) + ".fallback.part")
    n = 0
    with gzip.open(tmp, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLS)
        for line in proc.stdout:
            values = line.rstrip("\n").split("\t")
            if len(values) != len(COLS):
                continue
            if values[6].split(".", 1)[0] != gene_id:
                continue
            writer.writerow(values)
            n += 1
    stderr = proc.stderr.read() if proc.stderr is not None else ""
    code = proc.wait()
    if code != 0:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"tabix fallback failed for {gene} {context}: {stderr}")
    tmp.replace(target)
    done.write_text(f"complete\t{n}\n")
    return target, n, url


def main() -> None:
    with MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    failures = [r for r in rows if r["status"] == "TABIX_FAILED"]
    print(f"Recovering {len(failures)} failed GTEx requests", flush=True)
    recovered: dict[tuple[str, str], tuple[Path, int, str]] = {}
    for i, row in enumerate(failures, 1):
        print(f"[{i}/{len(failures)}] {row['gene']} | {row['context']}", flush=True)
        recovered[(row["gene"], row["context"])] = fetch(row)

    for row in rows:
        key = (row["gene"], row["context"])
        if key not in recovered:
            continue
        path, n, url = recovered[key]
        row["source"] = "GTEx_V8_imported_eQTL_Catalogue_fallback"
        row["source_url"] = url
        row["output"] = str(path)
        row["n_qtl_rows"] = str(n)
        row["status"] = "RECOVERED_FROM_OFFICIAL_GTEX_V8" if n else "NO_TARGET_GENE_QTL"
        row["error"] = ""

    fields = list(rows[0])
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    unresolved = [r for r in rows if r["status"] == "TABIX_FAILED"]
    print(f"Updated manifest: {MANIFEST}; unresolved={len(unresolved)}")
    if unresolved:
        raise RuntimeError(f"{len(unresolved)} requests remain unresolved")


if __name__ == "__main__":
    main()
