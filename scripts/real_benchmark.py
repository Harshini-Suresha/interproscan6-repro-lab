#!/usr/bin/env python3
"""real_benchmark.py

Measure reproducibility and cost of the real annotation data produced by
fetch_real_data.py, and validate it against a second real source, the InterPro
API.

This does not run InterProScan. It times the local normalise-and-fingerprint
step (the same step a real run would use) and it checks the UniProt-derived
annotations against the InterPro API for the same proteins.

Usage: python3 scripts/real_benchmark.py
Writes: results/benchmark_runs.csv, results/reproducibility_summary.csv,
        results/reference_agreement.csv
"""
import hashlib
import json
import os
import random
import sys
import time
import tracemalloc
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
ANN = os.path.join(ROOT, "annotations")
RES = os.path.join(ROOT, "results")

DATASETS = [
    ("coli", "UP000000625", "Escherichia coli"),
    ("mtb", "UP000001584", "Mycobacterium tuberculosis"),
]
CONFIGS = [("chunk1000_in_order", 1000, "in"), ("chunk100_in_order", 100, "in"),
           ("chunk10_shuffled", 10, "shuffle"), ("chunk500_reversed", 500, "reverse")]
REPEATS = 3
SAMPLE = 100  # proteins compared against the InterPro API per dataset


def http_json(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "interproscan6-repro-lab/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read()
            return json.loads(body) if body.strip() else {}
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if i < tries - 1:
                time.sleep(2 * (i + 1))
                continue
            raise


def interpro_api_entries(acc):
    """Set of InterPro accessions reported by the InterPro API for one protein."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/?page_size=200"
    iprs = set()
    while url:
        data = http_json(url)
        for r in data.get("results", []):
            if r.get("metadata", {}).get("source_database") == "interpro":
                a = r["metadata"].get("accession")
                if a:
                    iprs.add(a)
        url = data.get("next")
    return iprs


def bare_accession(header_token):
    parts = header_token.split("|")
    return parts[1] if len(parts) == 3 else parts[0]


def parse_fasta(path):
    seqs, name = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs.append(name)
                name = bare_accession(line[1:].strip().split()[0])
    if name:
        seqs.append(name)
    return seqs


def fingerprint(lines):
    """Drop the date column (11), sort, sha256; first 12 hex chars."""
    norm = []
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 15:
            norm.append("\t".join(f[:10] + f[11:15]))
        else:
            norm.append(ln.rstrip("\n"))
    norm.sort()
    return hashlib.sha256("\n".join(norm).encode()).hexdigest()[:12]


def peak_mb():
    return tracemalloc.get_traced_memory()[1] / (1024 * 1024)


def run_config(lines, chunk, order):
    work = list(lines)
    if order == "reverse":
        work.reverse()
    elif order == "shuffle":
        random.Random(42).shuffle(work)
    tracemalloc.reset_peak()
    t0 = time.perf_counter()
    seen = []
    for i in range(0, len(work), chunk):
        seen.extend(work[i:i + chunk])
    fp = fingerprint(seen)
    dt = (time.perf_counter() - t0) * 1000
    return dt, peak_mb(), fp


def annotated_proteins(path):
    prots = set()
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 5 and f[0]:
            prots.add(f[0])
    return prots


def protein_ipr_pairs(path):
    """(protein, InterPro accession) pairs from column 12."""
    pairs = set()
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 12 and f[0] and f[11] not in ("-", ""):
            pairs.add((f[0], f[11]))
    return pairs


def main():
    tracemalloc.start()
    bench, repro, cmprows = [], [], []
    for slug, proteome, species in DATASETS:
        ann = os.path.join(ANN, f"{slug}_{proteome}_annotations.tsv")
        fasta = os.path.join(DATA, f"{slug}_{proteome}.fasta")
        lines = open(ann).read().splitlines()
        accs = parse_fasta(fasta)
        n = len(accs)
        n_ann = len(annotated_proteins(ann))
        fps = []
        for cfg, chunk, order in CONFIGS:
            for rep in range(1, REPEATS + 1):
                dt, mem, fp = run_config(lines, chunk, order)
                fps.append(fp)
                bench.append({
                    "dataset": slug, "n_proteins": n, "source": f"UniProt {proteome}",
                    "config": cfg, "repeat": rep, "wall_ms": round(dt, 1),
                    "peak_heap_mb": round(mem, 1), "matches": len(lines),
                    "coverage_pct": round(100.0 * n_ann / n, 2),
                    "fingerprint": fp,
                })
        repro.append({
            "dataset": slug, "runs": len(fps), "distinct_fingerprints": len(set(fps)),
            "reproducible": str(len(set(fps)) == 1).lower(),
        })

        # cross-source agreement: UniProt-derived annotations vs InterPro API
        sample = sorted(accs)[:SAMPLE]
        local = protein_ipr_pairs(ann)
        uni_pairs = {(p, i) for p, i in local if p in set(sample)}
        api_pairs = set()
        for p in sample:
            for i in interpro_api_entries(p):
                api_pairs.add((p, i))
        shared = uni_pairs & api_pairs
        union = uni_pairs | api_pairs
        cmprows.append({
            "dataset": slug, "n_sampled": len(sample),
            "uniprot_derived_pairs": len(uni_pairs), "interpro_api_pairs": len(api_pairs),
            "shared": len(shared),
            "agreement_pct": round(100.0 * len(shared) / len(union), 1) if union else 0.0,
        })
        sys.stderr.write(f"{slug}: {n} proteins, {n_ann} annotated; "
                         f"cross-source {len(shared)}/{len(union)} pairs agree\n")

    _write(os.path.join(RES, "benchmark_runs.csv"), bench)
    _write(os.path.join(RES, "reproducibility_summary.csv"), repro)
    _write(os.path.join(RES, "reference_agreement.csv"), cmprows)


def _write(path, rows):
    cols = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(_q(r[c]) for c in cols) + "\n")


def _q(v):
    s = str(v)
    return '"' + s.replace('"', '""') + '"' if ("," in s or '"' in s) else s


if __name__ == "__main__":
    main()
