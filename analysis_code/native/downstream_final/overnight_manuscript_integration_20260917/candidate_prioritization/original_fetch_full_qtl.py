#!/usr/bin/env python3
"""Fetch complete +/-1-Mb cis-eQTL statistics for all retained Tier 1 genes.

The primary panel comprises seven vascular-relevant GTEx v8 contexts for every
gene. Four genes nominated in macrophage states additionally receive four
Alasoo_2018 macrophage contexts and two BLUEPRINT myeloid contexts. Remote tabix
requests are deliberately serial and rate-limited as recommended by EMBL-EBI.
An empty/low-expression gene-context is retained in the manifest rather than
silently discarded.
"""

from __future__ import annotations

import csv
import concurrent.futures
import gzip
import subprocess
import time
from pathlib import Path


ROOT = Path("/Users/taki/Desktop/Codex任务1/Tier1_full_colocalisation_v1")
COORDS = ROOT / "00_metadata" / "tier1_gene_coordinates.tsv"
OUT = ROOT / "01_qtl_regions_uniform"
MANIFEST = ROOT / "00_metadata" / "qtl_download_manifest.tsv"
TABIX = Path("/opt/homebrew/bin/tabix")
RATE_LIMIT_SECONDS = 2.25
REQUEST_TIMEOUT_SECONDS = 360
MAX_ATTEMPTS = 3

SUMSTATS_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"

GTEX = {
    "Artery_Aorta": ("QTD000131", 387),
    "Artery_Coronary": ("QTD000136", 213),
    "Artery_Tibial": ("QTD000141", 584),
    "Cells_Cultured_fibroblasts": ("QTD000216", 483),
    "Heart_Atrial_Appendage": ("QTD000251", 372),
    "Heart_Left_Ventricle": ("QTD000256", 382),
    "Whole_Blood": ("QTD000356", 670),
}

IMMUNE = {
    "Alasoo_macrophage_naive": ("QTS000001", "QTD000001", 84, "GRCh38"),
    "Alasoo_macrophage_IFNg": ("QTS000001", "QTD000006", 84, "GRCh38"),
    "Alasoo_macrophage_Salmonella": ("QTS000001", "QTD000011", 84, "GRCh38"),
    "Alasoo_macrophage_IFNg_Salmonella": ("QTS000001", "QTD000016", 84, "GRCh38"),
    "BLUEPRINT_monocyte": ("QTS000002", "QTD000021", 191, "GRCh37"),
    "BLUEPRINT_neutrophil": ("QTS000002", "QTD000026", 196, "GRCh37"),
}

GTEX_COLS = [
    "variant", "r2", "pvalue", "molecular_trait_object_id",
    "molecular_trait_id", "maf", "gene_id", "median_tpm", "beta", "se",
    "an", "ac", "chromosome", "position", "ref", "alt", "type", "rsid",
]

CATALOGUE_COLS = [
    "molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
    "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid",
]


def read_genes() -> list[dict[str, str]]:
    with COORDS.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def output_path(gene: str, context: str) -> Path:
    return OUT / f"{gene}__{context}__cis1Mb.tsv.gz"


def fetch_one(gene: dict[str, str], context: str, source: str,
              url: str, sample_size: int, build: str) -> dict[str, object]:
    path = output_path(gene["gene"], context)
    done = Path(str(path) + ".done")
    if path.exists() and done.exists() and path.stat().st_size > 70:
        done_fields = done.read_text().strip().split("\t")
        n_rows = int(done_fields[1]) if len(done_fields) > 1 else None
        return {
            "gene": gene["gene"], "gene_id": gene["gene_id"],
            "nominated_lineage": gene["nominated_lineage"], "context": context,
            "source": source, "sample_size": sample_size, "source_build": build,
            "source_url": url, "output": str(path), "n_qtl_rows": n_rows,
            "status": "REUSED" if n_rows else "NO_TARGET_GENE_QTL",
        }

    start_key = "start_grch37" if build == "GRCh37" else "start_grch38"
    end_key = "end_grch37" if build == "GRCh37" else "end_grch38"
    if gene[start_key] in {"", "NA"} or gene[end_key] in {"", "NA"}:
        return {
            "gene": gene["gene"], "gene_id": gene["gene_id"],
            "nominated_lineage": gene["nominated_lineage"], "context": context,
            "source": source, "sample_size": sample_size, "source_build": build,
            "source_url": url, "output": "", "n_qtl_rows": 0,
            "status": "SKIPPED_NO_BUILD_COORDINATE",
        }

    chrom = gene["chromosome"]
    start = max(1, int(gene[start_key]) - 1_000_000)
    end = int(gene[end_key]) + 1_000_000
    region = f"{chrom}:{start}-{end}"
    tmp = Path(str(path) + ".part")
    columns = CATALOGUE_COLS

    command = [str(TABIX), url, region]
    completed = None
    last_error = ""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            completed = subprocess.run(
                command, capture_output=True, text=True,
                timeout=REQUEST_TIMEOUT_SECONDS, check=False,
            )
            if completed.returncode == 0:
                break
            last_error = completed.stderr.strip()
        except subprocess.TimeoutExpired:
            last_error = f"timeout after {REQUEST_TIMEOUT_SECONDS}s"
        if attempt < MAX_ATTEMPTS:
            time.sleep(5 * attempt)

    if completed is None or completed.returncode != 0:
        tmp.unlink(missing_ok=True)
        return {
            "gene": gene["gene"], "gene_id": gene["gene_id"],
            "nominated_lineage": gene["nominated_lineage"], "context": context,
            "source": source, "sample_size": sample_size, "source_build": build,
            "source_url": url, "region_query": region, "output": "",
            "n_qtl_rows": 0, "status": "TABIX_FAILED", "error": last_error,
        }

    kept = 0
    with gzip.open(tmp, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GTEX_COLS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for line in completed.stdout.splitlines():
            values = line.rstrip("\n").split("\t")
            if len(values) != len(columns):
                continue
            record = dict(zip(columns, values))
            if record.get("gene_id", "").split(".", 1)[0] != gene["gene_id"]:
                continue
            writer.writerow({key: record.get(key, "") for key in GTEX_COLS})
            kept += 1
    tmp.replace(path)
    done.write_text(f"complete\t{kept}\n")
    return {
        "gene": gene["gene"], "gene_id": gene["gene_id"],
        "nominated_lineage": gene["nominated_lineage"], "context": context,
        "source": source, "sample_size": sample_size, "source_build": build,
        "source_url": url, "region_query": region, "output": str(path),
        "n_qtl_rows": kept,
        "status": "DOWNLOADED" if kept else "NO_TARGET_GENE_QTL",
    }


def main() -> None:
    if not TABIX.exists():
        raise FileNotFoundError(TABIX)
    OUT.mkdir(parents=True, exist_ok=True)
    genes = read_genes()
    jobs: list[tuple[dict[str, str], str, str, str, int, str]] = []
    for gene in genes:
        for context, (dataset, n) in GTEX.items():
            url = f"{SUMSTATS_BASE}/QTS000015/{dataset}/{dataset}.all.tsv.gz"
            jobs.append((gene, context, "GTEx_eQTL_Catalogue", url, n, "GRCh38"))
        if gene["nominated_lineage"] == "Macrophage":
            for context, (study, dataset, n, build) in IMMUNE.items():
                url = f"{SUMSTATS_BASE}/{study}/{dataset}/{dataset}.all.tsv.gz"
                jobs.append((gene, context, study.replace("QTS000001", "Alasoo_2018").replace("QTS000002", "BLUEPRINT"), url, n, build))

    rows: list[dict[str, object]] = []
    pending: dict[concurrent.futures.Future, tuple[int, tuple]] = {}
    next_job = 0
    last_start = 0.0
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        while next_job < len(jobs) or pending:
            while next_job < len(jobs) and len(pending) < 4:
                wait = RATE_LIMIT_SECONDS - (time.monotonic() - last_start)
                if wait > 0:
                    time.sleep(wait)
                job = jobs[next_job]
                gene, context, *_ = job
                print(f"[{next_job + 1}/{len(jobs)}] {gene['gene']} | {context}", flush=True)
                future = pool.submit(fetch_one, *job)
                pending[future] = (next_job + 1, job)
                next_job += 1
                last_start = time.monotonic()
            done, _ = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                _, job = pending.pop(future)
                row = future.result()
                rows.append(row)
                print(
                    f"  DONE {row['gene']} | {row['context']} | "
                    f"{row['status']} n={row.get('n_qtl_rows', 0)}",
                    flush=True,
                )

    fields = sorted({key for row in rows for key in row})
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    failed = [r for r in rows if r["status"] == "TABIX_FAILED"]
    print(f"Manifest: {MANIFEST}")
    print(f"Jobs={len(rows)} failed={len(failed)}")
    if failed:
        raise RuntimeError(f"{len(failed)} remote tabix requests failed")


if __name__ == "__main__":
    main()
