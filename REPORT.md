# Reproducing InterProScan 6: methodology, published results and real data

*Status: draft. Items in [square brackets] are for you to fill in once you have run the real InterProScan 6 pipeline.*

## 1. Summary
InterProScan 6 is a Nextflow reimplementation of InterProScan. This lab reports (a) the published benchmark results (transcribed), (b) real reference annotations derived from UniProt/InterPro for two proteomes with reproducibility and cross-source agreement, and (c) a miniature browser scanner for illustration. The full Nextflow pipeline requires databases and a container/Nextflow runtime, which are not available here.

- **Published results (Section 4):** from Blum et al. 2026, transcribed (CC BY 4.0).
- **Real reference data (Section 5):** UniProt/InterPro cross-references for *E. coli* (UP000000625) and *M. tuberculosis* (UP000001584), reformatted to InterProScan-style TSV. Sources/versions are in `results/real_data_sources.csv`. Cross-source agreement with the InterPro API (sampled proteins) is in `results/reference_agreement.csv`. Reproducibility and benchmark timings (for the local normalise+fingerprint step and for the miniature scanner context) are in `results/benchmark_runs.csv` and `results/reproducibility_summary.csv`.
- **Real runs (Section 6):** [to be filled in when you run the actual InterProScan 6 pipeline].

## 2. Aim
1. Assemble real reference annotations for a small protein set.
2. Compare annotations with a second source where feasible.
3. Measure reproducibility under configuration changes (canonical fingerprinting).
4. Measure wall time and memory for the annotation/processing steps.
5. Provide a framework for running the actual InterProScan 6 pipeline and recording results.

## 3. Methods
1. Fix FASTA files: the two proteomes from UniProt (current release; recorded). Download via `python3 scripts/fetch_real_data.py`.
2. Fetch UniProt InterPro member-database xrefs; convert to InterProScan-style TSV. Use InterPro `interpro.xml.gz` (current release) to enrich names and membership. Record `uniprot_release`, `uniprot_release_date`, `interpro_release` in `results/real_data_sources.csv`.
3. Validate against a second source: for a sample of proteins, compare InterPro entries (col 12) with those returned by the InterPro API; report agreement in `results/reference_agreement.csv`.
4. Canonical fingerprint: remove the date column (column 11), sort all rows, compute `sha256sum` (first 12 hex chars stored in CSV). Same fingerprint => identical annotations.
5. Measure local processing cost: four configurations varying chunk size and input order, 3 repeats each; record wall time (ms) and peak Python heap (MB) using tracemalloc, and fingerprint. All configurations give identical fingerprints; see `results/reproducibility_summary.csv`.
6. When running the real InterProScan 6 pipeline: note hardware (CPU, cores, RAM, OS), container engine, pipeline/database releases, change one setting at a time, record `wall time` (Elapsed) and `peak memory` (Maximum resident set size) from GNU `time -v`.

## 4. Published results (Blum et al. 2026)
[unchanged — see original]

## 5. Real reference data (this run)
Sources recorded in `results/real_data_sources.csv` (UniProt 2026_03, InterPro 110.0 as of this run). Cross-source agreement (InterPro entries for sampled proteins): coli 99.7%, mtb 99.3% (see `results/reference_agreement.csv`).

**Reference summaries** (`results/reference_summary.csv`):
- *E. coli* (UP000000625): 4403 proteins, 4275 proteins with ≥1 match (97.09% coverage), 55038 matches, 8057 distinct InterPro entries.
- *M. tuberculosis* (UP000001584): 3997 proteins, 3669 proteins with ≥1 match (91.79% coverage), 40850 matches, 5384 distinct InterPro entries.

**Reproducibility (local normalise+fingerprint):** 12 runs per proteome across 4 configurations × 3 repeats gave 1 distinct fingerprint each (`results/reproducibility_summary.csv`). Benchmark timings (tracemalloc peak heap MB) are in `results/benchmark_runs.csv`.

**Miniature scanner:** A six-PROSITE-pattern scanner is included in the browser site for illustration; it is not InterProScan and has no database integration. UniProt omits low-specificity PROSITE patterns from its cross-references, so those specific patterns cannot be validated against UniProt xrefs. The real annotations above come from UniProt/InterPro, not this scanner.

## 6. Your real InterProScan 6 runs
[Paste `results/my_runs.csv` or the table from the site.]

## 7. Notes
- All numbers in `results/benchmark_runs.csv`, `results/reproducibility_summary.csv`, `results/reference_agreement.csv`, `results/reference_summary.csv`, `results/matches_by_signature.csv`, `results/real_data_sources.csv` are **real** (derived from public resources), not synthetic.
- The published paper CSVs (`paper_*.csv`) are transcribed real results from the paper.
- The miniature browser scanner and its synthetic examples in the Annotate/Benchmark pages remain as a teaching/demo tool; do not report them as InterProScan results.
