#!/usr/bin/env python3
"""run_matrix.py (needs Python 3.9+, Linux or macOS)
run_matrix.py FASTA [--configs configs.txt] [--repeats 3] [--datadir data] [--out results/my_runs.csv]
configs.txt has one configuration per line: name|extra nextflow arguments
  docker|-profile docker
  singularity|-profile singularity
Check the InterProScan 6 docs for CPU and batch-size options and add them to the extra arguments.
Writes a CSV that the site's Your results page can import."""
import argparse, csv, glob, hashlib, os, shlex, subprocess, sys, time
ap = argparse.ArgumentParser()
ap.add_argument("fasta"); ap.add_argument("--configs"); ap.add_argument("--repeats", type=int, default=3)
ap.add_argument("--datadir", default="data"); ap.add_argument("--out", default="results/my_runs.csv")
ap.add_argument("--version", default="6.0.1"); ap.add_argument("--interpro", default="108.0")
ap.add_argument("--nextflow", default="nextflow")
a = ap.parse_args()
fasta, datadir = os.path.abspath(a.fasta), os.path.abspath(a.datadir)
cfgs = [("docker", "-profile docker")]
if a.configs:
    cfgs = [tuple(l.strip().split("|", 1)) for l in open(a.configs) if "|" in l and not l.startswith("#")]
n_in = sum(1 for l in open(fasta) if l.startswith(">"))
def pairs(path):
    s = set()
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 5 and f[4] not in ("", "-"): s.add((f[0], f[4]))
    return s
def fingerprint(path):
    rows = []
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        rows.append("\t".join(f[:10] + f[11:15]))
    return hashlib.sha256("\n".join(sorted(rows)).encode()).hexdigest()[:12]
def run_timed(cmd, cwd):
    """Returns (return code, wall minutes, peak GB of the largest process, via the OS)."""
    t0 = time.time(); pr = subprocess.Popen(cmd, cwd=cwd)
    _, status, ru = os.wait4(pr.pid, 0); pr.returncode = os.waitstatus_to_exitcode(status)
    kb = ru.ru_maxrss / (1024 if sys.platform == "darwin" else 1)  # macOS reports bytes, Linux kilobytes
    return pr.returncode, (time.time() - t0) / 60, kb / 1048576
H = ["Run name", "Configuration", "Input proteins", "Wall time (min)", "Peak memory (GB)",
     "Proteins annotated", "Agreement with reference (%)", "Output fingerprint", "Notes"]
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
rows, base = [], None
for name, extra in cfgs:
    for r in range(1, a.repeats + 1):
        wd = os.path.abspath(f"runs/{name}-r{r}"); os.makedirs(wd, exist_ok=True)
        cmd = [a.nextflow, "run", "ebi-pf-team/interproscan6",
               "-r", a.version, "--input", fasta, "--datadir", datadir, "--interpro", a.interpro] + shlex.split(extra)
        print("running", name, "repeat", r, flush=True)
        rc, mins, gb = run_timed(cmd, wd)
        if rc != 0:
            rows.append([f"{name}-r{r}", extra, n_in, "", "", "", "", "", "run failed, see runs/ folder"]); continue
        tsv = sorted(glob.glob(os.path.join(wd, "*.tsv")))[0]
        p = pairs(tsv); base = base if base is not None else p
        u = len(p | base)
        rows.append([f"{name}-r{r}", extra, n_in, f"{mins:.2f}", f"{gb:.2f}", len({x for x, _ in p}),
                     f"{100 * len(p & base) / u:.1f}" if u else "", fingerprint(tsv),
                     "agreement is against the first run, not a reference; memory is the largest host process only, so container memory is not included with docker"])
with open(a.out, "w", newline="") as f:
    w = csv.writer(f, quoting=csv.QUOTE_ALL); w.writerow(H); w.writerows(rows)
print("wrote", a.out)
