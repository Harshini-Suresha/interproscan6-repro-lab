#!/usr/bin/env python3
"""Build site/viz-data.js: data for the in-browser visualisation gallery.

Every value is computed here from the real files in this repository:
  annotations/coli_UP000000625_annotations.tsv   real UniProt/InterPro cross-references
  annotations/mtb_UP000001584_annotations.tsv    (InterProScan-style TSV, 15 columns)
  data/*.fasta                                   real UniProt proteome sequences
  results/*.csv                                  real summaries, benchmarks and the
                                                 transcribed paper tables (Blum et al. 2026)
No network access and no synthetic values.

Run:  python3 scripts/build_visualizations.py
"""

import csv
import json
import os
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATASETS = [
    {
        "key": "coli",
        "species": "E. coli",
        "proteome": "UP000000625",
        "fasta": "data/coli_UP000000625.fasta",
        "tsv": "annotations/coli_UP000000625_annotations.tsv",
        "color": "#0b6e6e",
        "pale": "#9dc6c6",
    },
    {
        "key": "mtb",
        "species": "M. tuberculosis",
        "proteome": "UP000001584",
        "fasta": "data/mtb_UP000001584.fasta",
        "tsv": "annotations/mtb_UP000001584_annotations.tsv",
        "color": "#b4541a",
        "pale": "#e0b197",
    },
]

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA_MASS = {  # average residue masses, Da
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
KD = {  # Kyte-Doolittle hydropathy
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def path(*parts):
    return os.path.join(ROOT, *parts)


def read_csv(rel):
    with open(path(rel), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_fasta(rel):
    seqs, cur, buf = [], None, []
    with open(path(rel), encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs.append((cur, "".join(buf)))
                cur, buf = line[1:].split()[0], []
                if cur.count("|") >= 2:  # sp|A5A616|NAME -> A5A616
                    cur = cur.split("|")[1]
            elif cur is not None:
                buf.append(line.strip())
    if cur is not None:
        seqs.append((cur, "".join(buf)))
    return seqs


def read_annotations(rel):
    rows = []
    with open(path(rel), encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 13:
                rows.append(f)
    return rows


def bin_counts(values, edges, labels):
    """labels has len(edges)+1 entries: label 0 is "< edges[0]",
    label i is [edges[i-1], edges[i]), the last label is edges[-1]..inf."""
    out = [0] * len(labels)
    for v in values:
        out[bisect_right(edges, v)] += 1
    return out


def pct(part, whole):
    return round(100.0 * part / whole, 2) if whole else 0.0


def sample(seq, k):
    if len(seq) <= k:
        return list(seq)
    step = (len(seq) - 1) / (k - 1)
    return [seq[round(i * step)] for i in range(k)]


def load_organism(meta):
    seqs = read_fasta(meta["fasta"])
    lens = {sid: len(s) for sid, s in seqs}
    aa = Counter()
    mw, gravy, charge = [], [], []
    for _, s in seqs:
        std = [c for c in s if c in AA_MASS]
        aa.update(std)
        if std:
            mw.append(sum(AA_MASS[c] for c in std) / 1000.0 + 0.0180153)
            gravy.append(sum(KD[c] for c in std) / len(std))
            charge.append(sum(1 for c in std if c in "DEKR") / len(std))

    rows = read_annotations(meta["tsv"])
    matches = Counter()
    db_match, db_prot, db_ipr = Counter(), defaultdict(set), Counter()
    sig_match, ipr_match, ipr_prot = Counter(), Counter(), defaultdict(set)
    ipr_set = defaultdict(set)
    for f in rows:
        prot, db, sig, ipr = f[0], f[3], f[4], f[11]
        matches[prot] += 1
        if db != "InterPro":
            db_match[db] += 1
            db_prot[db].add(prot)
            sig_match[sig] += 1
            if ipr not in ("-", "NULL", ""):
                db_ipr[db] += 1
        if ipr not in ("-", "NULL", ""):
            ipr_match[ipr] += 1
            ipr_prot[ipr].add(prot)
            ipr_set[prot].add(ipr)

    return {
        "meta": meta,
        "seqs": seqs,
        "lens": lens,
        "n": len(seqs),
        "aa": aa,
        "mw": mw,
        "gravy": gravy,
        "charge": charge,
        "matches": matches,
        "db_match": db_match,
        "db_prot": db_prot,
        "db_ipr": db_ipr,
        "sig_match": sig_match,
        "ipr_match": ipr_match,
        "ipr_prot": ipr_prot,
        "ipr_set": ipr_set,
        "n_matches": sum(matches.values()),
        "n_matched": sum(1 for sid, _ in seqs if matches.get(sid, 0)),
    }


# ---------------------------------------------------------------------------
# chart spec helpers: the browser renderer understands these small dicts
# ---------------------------------------------------------------------------

def vbar(labels, values, **kw):
    return {"type": "vbar", "labels": labels,
            "series": [{"name": "", "values": values, "color": kw.pop("color", None)}], **kw}


def hbar(labels, values, **kw):
    return {"type": "hbar", "labels": labels,
            "series": [{"name": "", "values": values, "color": kw.pop("color", None)}], **kw}


def lollipop(labels, values, **kw):
    return {"type": "lollipop", "labels": labels,
            "series": [{"name": "", "values": values, "color": kw.pop("color", None)}], **kw}


def grouped(labels, series, **kw):
    return {"type": "grouped", "labels": labels, "series": series, **kw}


def stacked(labels, series, **kw):
    return {"type": "stacked", "labels": labels, "series": series, **kw}


def line(labels=None, x=None, series=None, **kw):
    d = {"type": "line", "series": series or [], **kw}
    if labels is not None:
        d["labels"] = labels
    if x is not None:
        d["xvals"] = x
    return d


def scatter(series, **kw):
    return {"type": "scatter", "series": series, **kw}


def donut(slices, **kw):
    return {"type": "donut", "slices": slices, **kw}


def heatmap(rows, cols, matrix, **kw):
    return {"type": "heatmap", "rows": rows, "cols": cols, "matrix": matrix, **kw}


def rngchart(labels, mins, maxs, means, **kw):
    return {"type": "range", "labels": labels, "mins": mins, "maxs": maxs, "means": means, **kw}


def fmt_n(x):
    return f"{x:,}"


def build():
    charts = []

    def add(title, note, spec, wide=False):
        spec["wide"] = wide or spec.get("wide", False)
        charts.append({"title": title, "note": note, "chart": spec})

    orgs = {d["key"]: load_organism(d) for d in DATASETS}
    coli, mtb = orgs["coli"], orgs["mtb"]

    # ------------------------------------------------------------------ A.
    # per-organism annotation views (real UniProt/InterPro cross-references)
    for org in (coli, mtb):
        m = org["meta"]
        src = m["tsv"]
        tag = m["species"]
        n = org["n"]

        lens = [org["lens"][sid] for sid, _ in org["seqs"]]
        edges = list(range(0, 1300, 100))
        labels = [f"{edges[i]}–{edges[i+1]-1}" for i in range(len(edges) - 1)] + ["1200+"]
        add(f"{tag}: protein length distribution",
            f"All {fmt_n(n)} real sequences in {m['fasta']}.",
            vbar(labels, bin_counts(lens, edges[1:], labels), color=m["color"], yUnit="proteins"))

        mpp = [org["matches"].get(sid, 0) for sid, _ in org["seqs"]]
        e2 = [1, 2, 3, 4, 5, 6, 8, 11, 15, 21]
        l2 = ["0", "1", "2", "3", "4", "5", "6–7", "8–10", "11–14", "15–20", "21+"]
        add(f"{tag}: matches per protein",
            f"{fmt_n(org['n_matches'])} matches over {fmt_n(n)} proteins ({src}).",
            vbar(l2, bin_counts(mpp, e2, l2), color=m["color"], yUnit="proteins"))

        dbs = [d for d, _ in org["db_match"].most_common(12)]
        add(f"{tag}: matches by member database",
            f"Match rows per database, top {len(dbs)} ({src}).",
            hbar(dbs, [org["db_match"][d] for d in dbs], color=m["color"], xUnit="matches"))

        add(f"{tag}: proteins covered by each member database",
            f"Distinct proteins with at least one match in the database ({src}).",
            hbar(dbs, [len(org["db_prot"][d]) for d in dbs], color=m["color"], xUnit="proteins"))

        total = sum(org["db_match"].values())
        slices = [{"name": d, "value": org["db_match"][d], "color": c}
                  for d, c in zip(dbs[:8], ["#0b6e6e", "#b4541a", "#3b4fa0", "#7a4fbf",
                                            "#0e8f6f", "#c23b73", "#8a6d1c", "#2b7fb8"])]
        slices.append({"name": "other", "value": total - sum(s["value"] for s in slices),
                       "color": "#c9d4dc"})
        add(f"{tag}: share of matches by database",
            f"All {fmt_n(total)} member-database matches, by database ({src}).",
            donut(slices, center=f"{fmt_n(total)}", centerUnit="matches", color=m["color"]))

        top = org["ipr_match"].most_common(15)
        add(f"{tag}: top InterPro entries by matches",
            f"Most-hit InterPro entries across {fmt_n(org['n_matches'])} matches ({src}).",
            hbar([ipr for ipr, _ in top],
                 [c for _, c in top], color=m["color"], xUnit="matches", wide=True))

        tsig = org["sig_match"].most_common(15)
        add(f"{tag}: top signatures by matches",
            f"Most-hit member-database signatures ({src}).",
            hbar([s for s, _ in tsig], [c for _, c in tsig],
                 color=m["color"], xUnit="matches", wide=True))

        freq = list(org["sig_match"].values())
        e3 = [2, 4, 8, 16, 32, 64, 128]
        l3 = ["1", "2–3", "4–7", "8–15", "16–31", "32–63", "64–127", "128+"]
        add(f"{tag}: how often each signature hits",
            f"{fmt_n(len(freq))} distinct signatures, binned by number of matches ({src}).",
            vbar(l3, bin_counts(freq, e3, l3), color=m["color"], yUnit="signatures", logY=True))

        ordered = sorted(freq, reverse=True)
        total_m = sum(ordered)
        xs = list(range(0, 101))
        ys, run = [], 0
        for xp in xs:
            idx = min(len(ordered) - 1, round(xp / 100 * (len(ordered) - 1)))
            run = sum(ordered[:idx + 1]) if xp == 0 else max(run, sum(ordered[:idx + 1]))
            ys.append(round(100.0 * sum(ordered[:idx + 1]) / total_m, 2))
        add(f"{tag}: match concentration across signatures",
            "Share of all matches contributed by the most-hit signatures (ranked order).",
            line(labels=[str(x) for x in xs], series=[
                {"name": "% of matches", "values": ys, "color": m["color"]}],
                xUnit="% of signatures", yUnit="% of matches", wide=True))

        iprp = [len(org["ipr_set"].get(sid, ())) for sid, _ in org["seqs"]]
        e4 = [1, 2, 3, 4, 5, 6, 8, 12]
        l4 = ["0", "1", "2", "3", "4", "5", "6–7", "8–11", "12+"]
        add(f"{tag}: distinct InterPro entries per protein",
            f"Counted over all {fmt_n(n)} proteins ({src}).",
            vbar(l4, bin_counts(iprp, e4, l4), color=m["color"], yUnit="proteins"))

        add(f"{tag}: share of matches integrated into InterPro",
            "Per member database: matches that carry an InterPro entry, % of its matches "
            f"({src}).",
            hbar(dbs, [pct(org["db_ipr"][d], org["db_match"][d]) for d in dbs],
                 color=m["color"], xUnit="%"))

        good, bad = org["n_matched"], n - org["n_matched"]
        add(f"{tag}: proteins with any annotation",
            f"{fmt_n(good)} of {fmt_n(n)} proteins have at least one match ({pct(good, n)}%).",
            donut([{"name": "≥1 match", "value": good, "color": m["color"]},
                   {"name": "no match", "value": bad, "color": "#d7e0e7"}],
                  center=pct(good, n), centerUnit="%", color=m["color"]))

        pts = [[org["lens"][sid], org["matches"].get(sid, 0)] for sid, _ in org["seqs"]]
        pts = sample(pts, 650)
        add(f"{tag}: protein length vs match count",
            f"Sampled {len(pts)} of {fmt_n(n)} real proteins ({src}); a point per protein.",
            scatter([{"name": tag, "points": pts, "color": m["color"]}],
                    xUnit="length (aa)", yUnit="matches", dots=True))

        be = [0, 100, 200, 300, 400, 500, 600, 800, 1000, 1400]
        bl = [f"{be[i]}–{be[i+1]}" for i in range(len(be) - 1)] + ["1400+"]
        sums = [0.0] * len(bl)
        cnts = [0] * len(bl)
        for sid, _ in org["seqs"]:
            L = org["lens"][sid]
            i = min(len(bl) - 1, max(0, bisect_left(be, L) - 1)) if L >= be[0] else 0
            sums[i] += org["matches"].get(sid, 0)
            cnts[i] += 1
        means = [round(sums[i] / cnts[i], 2) if cnts[i] else 0 for i in range(len(bl))]
        add(f"{tag}: mean matches per protein by length",
            f"Length bins of 100–400 aa over all {fmt_n(n)} proteins ({src}).",
            line(labels=bl, series=[{"name": "mean", "values": means,
                                     "color": m["color"]}], yUnit="matches"))

        prot_per_ipr = [len(v) for v in org["ipr_prot"].values()]
        e5 = [2, 4, 8, 16, 32, 64, 128, 256]
        l5 = ["1", "2–3", "4–7", "8–15", "16–31", "32–63", "64–127", "128–255", "256+"]
        add(f"{tag}: proteins per InterPro entry",
            f"{fmt_n(len(org['ipr_prot']))} distinct entries, binned by number of "
            f"proteins they annotate ({src}).",
            vbar(l5, bin_counts(prot_per_ipr, e5, l5), color=m["color"],
                 yUnit="entries", logY=True))

        rf = sorted(freq, reverse=True)[:1500]
        add(f"{tag}: signature rank–frequency",
            "Log–log plot: signatures ranked by match count (rank 1 left).",
            line(x=list(range(1, len(rf) + 1)),
                 series=[{"name": "matches", "values": rf, "color": m["color"]}],
                 xlog=True, ylog=True, dots=False,
                 xUnit="rank", yUnit="matches"))

    # ---------------------------------------------------------- A (cross)
    add("Coverage: any match vs InterPro entry",
        "Share of each real proteome with at least one member-database match or "
        "InterPro entry (results/reference_summary.csv).",
        grouped(["≥1 member-DB match", "≥1 InterPro entry"],
                [{"name": "E. coli", "values": [pct(coli["n_matched"], coli["n"]),
                                                pct(sum(1 for v in coli["ipr_set"].values()), coli["n"])],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis", "values": [pct(mtb["n_matched"], mtb["n"]),
                                                        pct(sum(1 for v in mtb["ipr_set"].values()), mtb["n"])],
                  "color": mtb["meta"]["color"]}],
                yUnit="%", yMax=100))

    e2b = [1, 2, 3, 4, 5, 6, 8, 11, 15, 21]
    l2b = ["0", "1", "2", "3", "4", "5", "6–7", "8–10", "11–14", "15–20", "21+"]
    a = [pct(c, coli["n"]) for c in bin_counts(
        [coli["matches"].get(sid, 0) for sid, _ in coli["seqs"]], e2b, l2b)]
    b = [pct(c, mtb["n"]) for c in bin_counts(
        [mtb["matches"].get(sid, 0) for sid, _ in mtb["seqs"]], e2b, l2b)]
    add("Matches per protein: E. coli vs M. tuberculosis",
        "Same bins, normalised to % of each proteome (real annotation TSVs).",
        grouped(l2b, [{"name": "E. coli", "values": a, "color": coli["meta"]["color"]},
                      {"name": "M. tuberculosis", "values": b, "color": mtb["meta"]["color"]}],
                yUnit="% of proteome", wide=True))

    e1b = list(range(0, 1300, 100))
    l1b = [f"{e1b[i]}–{e1b[i+1]-1}" for i in range(len(e1b) - 1)] + ["1200+"]
    a = [pct(c, coli["n"]) for c in bin_counts(
        [coli["lens"][sid] for sid, _ in coli["seqs"]], e1b[1:], l1b)]
    b = [pct(c, mtb["n"]) for c in bin_counts(
        [mtb["lens"][sid] for sid, _ in mtb["seqs"]], e1b[1:], l1b)]
    add("Length distribution: E. coli vs M. tuberculosis",
        "Real FASTA sequences, normalised to % of each proteome.",
        grouped(l1b, [{"name": "E. coli", "values": a, "color": coli["meta"]["color"]},
                      {"name": "M. tuberculosis", "values": b, "color": mtb["meta"]["color"]}],
                yUnit="% of proteome", wide=True))

    all_db = sorted(set(coli["db_match"]) | set(mtb["db_match"]),
                    key=lambda d: -(coli["db_match"][d] + mtb["db_match"][d]))
    add("Member-database matches: the two proteomes",
        "Real match counts per database in each annotation file.",
        grouped(all_db, [{"name": "E. coli", "values": [coli["db_match"][d] for d in all_db],
                          "color": coli["meta"]["color"]},
                         {"name": "M. tuberculosis",
                          "values": [mtb["db_match"][d] for d in all_db],
                          "color": mtb["meta"]["color"]}],
                yUnit="matches", wide=True))

    add("Member-database coverage of each proteome",
        "% of proteome proteins with ≥1 match in the database (real TSVs).",
        grouped(all_db, [{"name": "E. coli",
                          "values": [pct(len(coli["db_prot"][d]), coli["n"]) for d in all_db],
                          "color": coli["meta"]["color"]},
                         {"name": "M. tuberculosis",
                          "values": [pct(len(mtb["db_prot"][d]), mtb["n"]) for d in all_db],
                          "color": mtb["meta"]["color"]}],
                yUnit="% of proteome", wide=True))

    add("Total real matches in this project",
        "55,038 + 40,850 member-database-and-InterPro matches across both annotation files.",
        stacked(["all matches"],
                [{"name": "E. coli", "values": [coli["n_matches"]], "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis", "values": [mtb["n_matches"]], "color": mtb["meta"]["color"]}],
                yUnit="matches"))

    # ------------------------------------------------------------------ B.
    # results/matches_by_signature.csv (real signature-level frequencies)
    sig = read_csv("results/matches_by_signature.csv")
    by_ds = {"coli": [r for r in sig if r["dataset"] == "coli"],
             "mtb": [r for r in sig if r["dataset"] == "mtb"]}

    top = sorted(sig, key=lambda r: -int(r["matches"]))[:20]
    add("Top signatures across both proteomes",
        "Matches per signature from results/matches_by_signature.csv (top 20).",
        hbar([r["signature"] for r in top], [int(r["matches"]) for r in top],
             xUnit="matches", wide=True))

    db_names = sorted({r["member_database"] for r in sig},
                      key=lambda d: -sum(1 for r in sig if r["member_database"] == d))
    add("Distinct signatures per member database",
        f"{fmt_n(len(sig))} signature rows in results/matches_by_signature.csv.",
        vbar(db_names, [sum(1 for r in sig if r["member_database"] == d) for d in db_names],
             yUnit="signatures", wide=True))

    freq = [int(r["matches"]) for r in sig]
    e = [2, 4, 8, 16, 32, 64, 128, 256]
    l = ["1", "2–3", "4–7", "8–15", "16–31", "32–63", "64–127", "128–255", "256+"]
    add("Matches per signature (both proteomes)",
        "Distribution over all real signature rows.",
        vbar(l, bin_counts(freq, e, l), yUnit="signatures", logY=True))

    pfreq = [int(r["proteins_with_hit"]) for r in sig]
    add("Proteins per signature (both proteomes)",
        "How many proteins each real signature hits.",
        vbar(l, bin_counts(pfreq, e, l), yUnit="signatures", logY=True))

    for key in ("coli", "mtb"):
        rows = by_ds[key]
        m = coli if key == "coli" else mtb
        pts = [[int(r["proteins_with_hit"]), int(r["matches"])] for r in rows]
        add(f"{m['meta']['species']}: signature matches vs proteins hit",
            f"One point per real signature ({fmt_n(len(pts))} points); diagonal = every "
            "match in a different protein.",
            scatter([{"name": m["meta"]["species"], "points": sample(pts, 900),
                      "color": m["meta"]["color"]}],
                    xUnit="proteins with hit", yUnit="matches", diag=True))

        cand = [r for r in rows if int(r["proteins_with_hit"]) >= 10]
        cand.sort(key=lambda r: -(int(r["matches"]) / max(1, int(r["proteins_with_hit"]))))
        t = cand[:12]
        add(f"{m['meta']['species']}: signatures with most matches per protein",
            "Multi-domain repeats: signatures hitting the same protein several times "
            "(≥10 proteins).",
            hbar([r["signature"] for r in t],
                 [round(int(r["matches"]) / max(1, int(r["proteins_with_hit"])), 2) for r in t],
                 color=m["meta"]["color"], xUnit="matches per protein"))

    rate = []
    for d in db_names:
        rows = [r for r in sig if r["member_database"] == d]
        mapped = sum(1 for r in rows if r["interpro"] not in ("-", "", "NULL"))
        rate.append(pct(mapped, len(rows)))
    add("Signatures mapped to InterPro, by database",
        "Share of each database's signatures that carry an InterPro entry.",
        hbar(db_names, rate, xUnit="%"))

    ipr_sigs = Counter()
    for r in sig:
        if r["interpro"] not in ("-", "", "NULL"):
            ipr_sigs[r["interpro"]] += 1
    t = ipr_sigs.most_common(15)
    add("InterPro entries with the most signatures",
        "Distinct member-database signatures grouped under each InterPro entry.",
        hbar([k for k, _ in t], [v for _, v in t], xUnit="signatures", wide=True))

    sc = {r["signature"] for r in by_ds["coli"]}
    sm = {r["signature"] for r in by_ds["mtb"]}
    shared = len(sc & sm)
    add("Signatures shared between the two proteomes",
        "Signature IDs appearing in both real annotation sets.",
        donut([{"name": "both", "value": shared, "color": "#0b6e6e"},
               {"name": "E. coli only", "value": len(sc - sm), "color": "#b4541a"},
               {"name": "M. tuberculosis only", "value": len(sm - sc), "color": "#3b4fa0"}],
              center=fmt_n(shared), centerUnit="shared"))

    add("Signatures per database: the two proteomes",
        "Distinct signature IDs per database in each real annotation set.",
        grouped(db_names,
                [{"name": "E. coli",
                  "values": [sum(1 for r in by_ds["coli"] if r["member_database"] == d) for d in db_names],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis",
                  "values": [sum(1 for r in by_ds["mtb"] if r["member_database"] == d) for d in db_names],
                  "color": mtb["meta"]["color"]}],
                yUnit="signatures", wide=True))

    # ------------------------------------------------------------------ C.
    # paper Table 1: runtime, nine real reference proteomes (transcribed)
    t1 = read_csv("results/paper_table1_runtimes.csv")
    names = [r["proteome"] for r in t1]
    seqs = [int(r["sequences"]) for r in t1]
    i5 = [int(r["ips5_min"]) for r in t1]
    i6 = [int(r["ips6_min"]) for r in t1]
    lk = [int(r["ips6_plus_match_lookup_min"]) for r in t1]
    red = [float(r["ips6_reduction_pct_computed"]) for r in t1]
    spd = [float(r["speedup_x_computed"]) for r in t1]
    src = "results/paper_table1_runtimes.csv (Blum et al. 2026, Table 1)"

    add("Sequences per reference proteome", f"Nine UniProt reference proteomes; {src}.",
        vbar(names, seqs, yUnit="sequences", color="#3b4fa0", wide=True))
    add("InterProScan 5 runtime per proteome", f"Wall minutes, 8 CPUs; {src}.",
        vbar(names, i5, yUnit="minutes", color="#b4541a", wide=True))
    add("InterProScan 6 runtime per proteome", f"Wall minutes, 8 CPUs; {src}.",
        vbar(names, i6, yUnit="minutes", color="#0b6e6e", wide=True))
    add("Runtime with the match lookup switched on", f"Wall minutes; {src}.",
        vbar(names, lk, yUnit="minutes", color="#0e8f6f", wide=True))
    add("Minutes saved by InterProScan 6", "IPS 5 minus IPS 6 wall time (computed).",
        vbar(names, [round(a - b, 1) for a, b in zip(i5, i6)], yUnit="minutes",
             color="#0e8f6f", wide=True))
    add("Time saved by InterProScan 6, %", f"Computed reduction; {src}.",
        vbar(names, red, yUnit="%", color="#0b6e6e", wide=True))
    add("Speed-up from InterProScan 5 to 6", f"Ratio of wall times; {src}.",
        vbar(names, spd, yUnit="×", color="#7a4fbf", wide=True))
    add("IPS 5 vs IPS 6 side by side",
        f"Wall minutes for the same nine proteomes; {src}.",
        grouped(names, [{"name": "IPS 5", "values": i5, "color": "#b4541a"},
                        {"name": "IPS 6", "values": i6, "color": "#0b6e6e"}],
                yUnit="minutes", wide=True))
    add("IPS 6 vs match lookup",
        "Wall minutes with and without checksum lookup; note the scale change.",
        grouped(names, [{"name": "IPS 6", "values": i6, "color": "#0b6e6e"},
                        {"name": "with lookup", "values": lk, "color": "#8a6d1c"}],
                yUnit="minutes", wide=True))
    add("Milliseconds per sequence: IPS 5 vs IPS 6",
        "Wall time divided by proteome size (computed from Table 1).",
        grouped(names, [{"name": "IPS 5", "values": [round(a * 60000 / s, 1) for a, s in zip(i5, seqs)],
                         "color": "#b4541a"},
                        {"name": "IPS 6", "values": [round(a * 60000 / s, 1) for a, s in zip(i6, seqs)],
                         "color": "#0b6e6e"}],
                yUnit="ms per sequence", wide=True))

    for label, vals in (("InterProScan 5", i5), ("InterProScan 6", i6), ("match lookup", lk)):
        add(f"Proteome size vs {label} time",
            "Each point is one proteome; time grows with size (log–log).",
            scatter([{"name": label, "points": [[s, v] for s, v in zip(seqs, vals)],
                      "color": "#3b4fa0" if label.endswith("5") else "#0b6e6e"}],
                    xlog=True, ylog=True, xUnit="sequences", yUnit="minutes", dots=True))

    add("Reduction vs proteome size (bubble = speed-up)",
        "X: sequences (log). Y: % time saved. Bubble size: speed-up ×.",
        scatter([{"name": "proteome", "points": [[s, r, 5 + 3 * x] for s, r, x in zip(seqs, red, spd)],
                  "color": "#7a4fbf"}],
                xlog=True, xUnit="sequences (log)", yUnit="% saved"))

    order = sorted(range(len(names)), key=lambda i: red[i])
    add("Proteomes ranked by time saved",
        "Dot plot, largest reduction on the right.",
        lollipop([names[i] for i in reversed(order)], [red[i] for i in reversed(order)],
                 xUnit="% saved", color="#0b6e6e", wide=True))

    add("Total wall time across all nine proteomes",
        "Sum of Table 1 columns (computed): IPS 5 vs IPS 6 vs match lookup.",
        vbar(["IPS 5", "IPS 6", "with lookup"], [sum(i5), sum(i6), sum(lk)],
             yUnit="minutes", color="#3b4fa0"))

    # ------------------------------------------------------------------ D.
    # paper Table 2: agreement with InterProScan 5 (transcribed)
    t2 = read_csv("results/paper_table2_concordance.csv")
    dbs = [r["member_database"] for r in t2]
    prec = [float(r["precision"]) for r in t2]
    sens = [float(r["sensitivity"]) for r in t2]
    f1 = [float(r["f1"]) for r in t2]
    src2 = "results/paper_table2_concordance.csv (Blum et al. 2026, Table 2)"
    add("Precision per member database", f"Against IPS 5 on Swiss-Prot; {src2}.",
        lollipop(dbs, prec, xUnit="precision", color="#0b6e6e", wide=True))
    add("Sensitivity per member database", f"Against IPS 5 on Swiss-Prot; {src2}.",
        lollipop(dbs, sens, xUnit="sensitivity", color="#3b4fa0", wide=True))
    add("F1 per member database", f"Harmonic mean of precision and sensitivity; {src2}.",
        lollipop(dbs, f1, xUnit="F1", color="#7a4fbf", wide=True))
    add("Precision, sensitivity and F1 together", f"All thirteen databases; {src2}.",
        grouped(dbs, [{"name": "precision", "values": prec, "color": "#0b6e6e"},
                      {"name": "sensitivity", "values": sens, "color": "#3b4fa0"},
                      {"name": "F1", "values": f1, "color": "#b4541a"}],
                yUnit="score", yMin=0.9, wide=True))
    add("Precision vs sensitivity",
        "Zoomed axes; almost every database sits at the top-right corner (1, 1).",
        scatter([{"name": "database", "points": [[p, s] for p, s in zip(prec, sens)],
                  "color": "#7a4fbf"}],
                xUnit="precision", yUnit="sensitivity", diag=True, dots=True))
    add("Distance from perfect agreement (F1 = 1)",
        "F1 gap ×1000; only four databases are not exactly 1.000.",
        hbar(dbs, [round((1 - v) * 1000, 1) for v in f1],
             xUnit="gap ×1000", wide=True))
    add("Overall agreement ranking",
        "Mean of precision, sensitivity and F1 per database.",
        hbar(dbs, [round((p + s + v) / 3, 4) for p, s, v in zip(prec, sens, f1)],
             xUnit="mean score", wide=True))

    # ------------------------------------------------------------------ E.
    cloud = read_csv("results/paper_cloud_costs_human_proteome.csv")
    clabels = [f"{r['provider']} · {r['mode']}" for r in cloud]
    add("Cloud cost for one human proteome",
        "Section 3.2 of the paper; four real price points (USD).",
        vbar(clabels, [float(r["cost_usd"]) for r in cloud], yUnit="USD",
             color="#2b7fb8", wide=True))
    add("Single-instance vs hybrid batch cost",
        "Same four runs, grouped by mode.",
        grouped(["single instance", "hybrid (Batch, spot)"],
                [{"name": "AWS", "values": [10.14, 12.43], "color": "#b4541a"},
                 {"name": "Google Cloud", "values": [8.13, 11.53], "color": "#2b7fb8"}],
                yUnit="USD"))

    # ------------------------------------------------------------------ F.
    ref = {r["dataset"]: r for r in read_csv("results/reference_summary.csv")}
    rc, rm = ref["coli"], ref["mtb"]
    add("Proteins and proteins with a match",
        "Real reference summaries for the two proteomes (results/reference_summary.csv).",
        grouped(["proteins", "with ≥1 match"],
                [{"name": "E. coli", "values": [int(rc["n_proteins"]), int(rc["proteins_with_match"])],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis",
                  "values": [int(rm["n_proteins"]), int(rm["proteins_with_match"])],
                  "color": mtb["meta"]["color"]}],
                yUnit="proteins"))
    add("Coverage of each proteome",
        "Proteins with at least one match, % (results/reference_summary.csv).",
        vbar(["E. coli", "M. tuberculosis"],
             [float(rc["coverage_pct"]), float(rm["coverage_pct"])],
             yUnit="%", yMax=100, color="#0b6e6e"))
    add("Total matches per proteome",
        "All match rows in each real annotation file.",
        vbar(["E. coli", "M. tuberculosis"], [int(rc["total_matches"]), int(rm["total_matches"])],
             yUnit="matches", color="#3b4fa0"))
    add("Matches per covered protein",
        "Total matches divided by the number of proteins that have any (computed).",
        vbar(["E. coli", "M. tuberculosis"],
             [round(int(rc["total_matches"]) / int(rc["proteins_with_match"]), 2),
              round(int(rm["total_matches"]) / int(rm["proteins_with_match"]), 2)],
             yUnit="matches", color="#7a4fbf"))
    add("Distinct signatures vs InterPro entries",
        "How many distinct IDs the annotations use.",
        grouped(["signatures", "InterPro entries"],
                [{"name": "E. coli", "values": [int(rc["distinct_signatures"]),
                                                int(rc["distinct_interpro_entries"])],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis", "values": [int(rm["distinct_signatures"]),
                                                        int(rm["distinct_interpro_entries"])],
                  "color": mtb["meta"]["color"]}],
                yUnit="distinct IDs"))
    add("Matches per distinct signature",
        "Match rows divided by distinct signatures (computed); higher = more reuse.",
        vbar(["E. coli", "M. tuberculosis"],
             [round(int(rc["total_matches"]) / int(rc["distinct_signatures"]), 2),
              round(int(rm["total_matches"]) / int(rm["distinct_signatures"]), 2)],
             yUnit="matches per signature", color="#0e8f6f"))

    # ------------------------------------------------------------------ G.
    agr = read_csv("results/reference_agreement.csv")
    add("Agreement with the InterPro API",
        "Shared (protein, entry) pairs over everything seen, % (results/reference_agreement.csv).",
        vbar([r["dataset"] for r in agr],
             [float(r["agreement_pct"]) for r in agr],
             yUnit="%", yMin=98, yMax=100, color="#0b6e6e"))
    shared_pairs = [int(r["shared"]) for r in agr]
    diff_pairs = [int(r["uniprot_derived_pairs"]) - int(r["shared"]) + int(r["interpro_api_pairs"]) - int(r["shared"])
                  for r in agr]
    add("Where the two sources disagree",
        "Pairs shared by UniProt-derived and API results vs pairs seen in only one source.",
        stacked([r["dataset"] for r in agr],
                [{"name": "shared by both", "values": shared_pairs, "color": "#0b6e6e"},
                 {"name": "in one source only", "values": diff_pairs, "color": "#b4541a"}],
                yUnit="(protein, entry) pairs"))

    # ------------------------------------------------------------------ H.
    rep = read_csv("results/reproducibility_summary.csv")
    add("Reproducibility: runs vs distinct fingerprints",
        "12 runs each (4 configurations × 3 repeats) — one fingerprint per proteome.",
        grouped([r["dataset"] for r in rep],
                [{"name": "runs", "values": [int(r["runs"]) for r in rep], "color": "#3b4fa0"},
                 {"name": "distinct fingerprints", "values": [int(r["distinct_fingerprints"]) for r in rep],
                  "color": "#0b6e6e"}],
                yUnit="count"))

    # ------------------------------------------------------------------ I.
    bench = read_csv("results/benchmark_runs.csv")
    cfg_names = {"chunk1000_in_order": "chunk 1000", "chunk100_in_order": "chunk 100",
                 "chunk10_shuffled": "chunk 10 shuffle", "chunk500_reversed": "chunk 500 reverse"}
    cfgs = list(cfg_names)
    srcb = "results/benchmark_runs.csv (real local normalise + fingerprint step)"

    for key in ("coli", "mtb"):
        rows = [r for r in bench if r["dataset"] == key]
        m = coli if key == "coli" else mtb
        labels = [f"{cfg_names[r['config']]} #{r['repeat']}" for r in rows]
        add(f"{m['meta']['species']}: wall time of all 12 runs",
            f"Milliseconds per real benchmark run, 3 repeats per configuration; {srcb}.",
            vbar(labels, [float(r["wall_ms"]) for r in rows],
                 yUnit="ms", color=m["meta"]["color"], wide=True))

    means_w = {c: {k: round(sum(float(r["wall_ms"]) for r in bench
                                if r["dataset"] == k and r["config"] == c) / 3, 1)
                   for k in ("coli", "mtb")} for c in cfgs}
    add("Mean wall time per configuration",
        f"3 repeats each; {srcb}.",
        grouped([cfg_names[c] for c in cfgs],
                [{"name": "E. coli", "values": [means_w[c]["coli"] for c in cfgs],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis", "values": [means_w[c]["mtb"] for c in cfgs],
                  "color": mtb["meta"]["color"]}],
                yUnit="ms", wide=True))

    mins, maxs, means = [], [], []
    labs = []
    for k in ("coli", "mtb"):
        m = coli if k == "coli" else mtb
        for c in cfgs:
            v = [float(r["wall_ms"]) for r in bench if r["dataset"] == k and r["config"] == c]
            labs.append(f"{m['meta']['key']} · {cfg_names[c]}")
            mins.append(round(min(v), 1))
            maxs.append(round(max(v), 1))
            means.append(round(sum(v) / len(v), 1))
    add("Run-to-run spread inside each configuration",
        "Min–max lines with the mean as a dot; 3 real repeats per point.",
        rngchart(labs, mins, maxs, means, yUnit="ms", wide=True))

    add("Repeat spread (max − min) per configuration",
        "Lower is steadier; computed from the 3 real repeats.",
        grouped([cfg_names[c] for c in cfgs],
                [{"name": "E. coli",
                  "values": [round(max(float(r["wall_ms"]) for r in bench if r["dataset"] == "coli" and r["config"] == c) -
                                      min(float(r["wall_ms"]) for r in bench if r["dataset"] == "coli" and r["config"] == c), 1)
                             for c in cfgs], "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis",
                  "values": [round(max(float(r["wall_ms"]) for r in bench if r["dataset"] == "mtb" and r["config"] == c) -
                                      min(float(r["wall_ms"]) for r in bench if r["dataset"] == "mtb" and r["config"] == c), 1)
                             for c in cfgs], "color": mtb["meta"]["color"]}],
                yUnit="ms"))

    add("Peak heap per configuration",
        "Peak heap in MB reported per real run (identical within each proteome).",
        grouped([cfg_names[c] for c in cfgs],
                [{"name": "E. coli",
                  "values": [round(sum(float(r["peak_heap_mb"]) for r in bench if r["dataset"] == "coli" and r["config"] == c) / 3, 1)
                             for c in cfgs], "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis",
                  "values": [round(sum(float(r["peak_heap_mb"]) for r in bench if r["dataset"] == "mtb" and r["config"] == c) / 3, 1)
                             for c in cfgs], "color": mtb["meta"]["color"]}],
                yUnit="MB"))

    rates = {c: {k: round(int(next(r for r in bench if r["dataset"] == k)["n_proteins"]) /
                         (means_w[c][k] / 1000)) for k in ("coli", "mtb")} for c in cfgs}
    add("Processing rate per configuration",
        "Proteins per second, computed from the real mean wall times.",
        grouped([cfg_names[c] for c in cfgs],
                [{"name": "E. coli", "values": [rates[c]["coli"] for c in cfgs],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis", "values": [rates[c]["mtb"] for c in cfgs],
                  "color": mtb["meta"]["color"]}],
                yUnit="proteins/s"))

    add("Wall time vs peak heap across all runs",
        "One point per real run; colour marks the proteome.",
        scatter([{"name": "E. coli",
                  "points": [[float(r["wall_ms"]), float(r["peak_heap_mb"])]
                             for r in bench if r["dataset"] == "coli"],
                  "color": coli["meta"]["color"]},
                 {"name": "M. tuberculosis",
                  "points": [[float(r["wall_ms"]), float(r["peak_heap_mb"])]
                             for r in bench if r["dataset"] == "mtb"],
                  "color": mtb["meta"]["color"]}],
                xUnit="wall time (ms)", yUnit="peak heap (MB)", dots=True))

    add("Match count is identical in every run",
        "All 12 runs per proteome produced exactly the same number of matches; the lines "
        "are flat by design (sorted output is byte-identical too).",
        line(labels=[f"{i+1}" for i in range(12)],
             series=[{"name": "E. coli", "values": [55038] * 12, "color": coli["meta"]["color"]},
                     {"name": "M. tuberculosis", "values": [40850] * 12,
                      "color": mtb["meta"]["color"]}],
             yMin=38000, yMax=58000, xUnit="run", yUnit="matches", wide=True))

    add("Coverage is identical in every run",
        "97.09% and 91.79% in all configurations and repeats.",
        lollipop(["E. coli", "M. tuberculosis"], [97.09, 91.79],
                 xUnit="% coverage", color="#0b6e6e", yMin=0, yMax=100))

    # ------------------------------------------------------------------ J.
    # real sequences: composition and biophysics
    caa = [pct(coli["aa"][a], sum(coli["aa"].values())) for a in AA20]
    maa = [pct(mtb["aa"][a], sum(mtb["aa"].values())) for a in AA20]
    add("Amino-acid composition: E. coli",
        f"Counted over all {fmt_n(sum(coli['aa'].values()))} real residues.",
        vbar(list(AA20), caa, yUnit="%", color=coli["meta"]["color"], wide=True))
    add("Amino-acid composition: M. tuberculosis",
        f"Counted over all {fmt_n(sum(mtb['aa'].values()))} real residues.",
        vbar(list(AA20), maa, yUnit="%", color=mtb["meta"]["color"], wide=True))
    add("Amino-acid composition: the two proteomes",
        "Share of each residue, both real proteomes side by side.",
        grouped(list(AA20), [{"name": "E. coli", "values": caa, "color": coli["meta"]["color"]},
                             {"name": "M. tuberculosis", "values": maa,
                              "color": mtb["meta"]["color"]}],
                yUnit="% of residues", wide=True))
    add("Composition differences (M. tuberculosis − E. coli)",
        "Percentage-point difference per residue; positive bars are more common in Mtb.",
        vbar(list(AA20), [round(b - a, 2) for a, b in zip(caa, maa)],
             yUnit="difference (pp)", color="#3b4fa0", wide=True))

    heat = heatmap(["E. coli", "M. tuberculosis"], list(AA20), [caa, maa],
                   color="#0b6e6e", unit="%")
    add("Composition heat map", "Same real composition, as a map.",
        heat, wide=True)

    mwedges = [0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 160]
    mwlabels = [f"{mwedges[i]}–{mwedges[i+1]}" for i in range(len(mwedges) - 1)] + ["160+"]
    for org in (coli, mtb):
        add(f"{org['meta']['species']}: predicted molecular weight",
            "Per real protein, from average residue masses (kDa).",
            vbar(mwlabels, bin_counts(org["mw"], mwedges[1:], mwlabels),
                 color=org["meta"]["color"], yUnit="proteins", wide=True))

    gedges = [-1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4]
    glabels = ["<−1.4", "−1.4–−1.2", "−1.2–−1.0", "−1.0–−0.8", "−0.8–−0.6",
               "−0.6–−0.4", "−0.4–−0.2", "−0.2–0.0", "0.0–0.2", "0.2–0.4", "0.4+"]
    for org in (coli, mtb):
        add(f"{org['meta']['species']}: average hydrophobicity (GRAVY)",
            "Kyte–Doolittle mean per real protein; most proteins are mildly hydrophilic.",
            vbar(glabels, bin_counts(org["gravy"], gedges, glabels),
                 color=org["meta"]["color"], yUnit="proteins", wide=True))

    cedges = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
    clabels2 = ["0–4%", "5–9%", "10–14%", "15–19%", "20–24%", "25–29%", "30–34%", "35–39%", "40%+"]
    for org in (coli, mtb):
        add(f"{org['meta']['species']}: charged-residue fraction",
            "D+E+K+R as a share of each real protein; Mtb is known to be arginine-rich.",
            vbar(clabels2, bin_counts(org["charge"], cedges[1:], clabels2),
                 color=org["meta"]["color"], yUnit="proteins", wide=True))

    xs = [round(x * 50) for x in range(0, 101)]
    cl = sorted(coli["lens"].values())
    ml = sorted(mtb["lens"].values())
    yc = [pct(bisect_left(cl, x), len(cl)) for x in xs]
    ym = [pct(bisect_left(ml, x), len(ml)) for x in xs]
    add("Protein length CDF: E. coli vs M. tuberculosis",
        "Cumulative share of proteins at or below each length, real FASTAs.",
        line(x=xs, series=[{"name": "E. coli", "values": yc, "color": coli["meta"]["color"]},
                           {"name": "M. tuberculosis", "values": ym,
                            "color": mtb["meta"]["color"]}],
             xUnit="length (aa)", yUnit="% of proteome", wide=True))

    longest = sorted(coli["seqs"], key=lambda kv: -len(kv[1]))[:10]
    add("Ten longest E. coli proteins",
        "Real sequences from the proteome FASTA.",
        hbar([sid for sid, _ in longest], [len(s) for _, s in longest],
             xUnit="length (aa)", color="#3b4fa0", wide=True))

    # ------------------------------------------------------------------ K.
    # position-resolved composition of the real sequences: where in the chain
    # is each residue over- or under-represented relative to the whole proteome
    NB = 10
    for org in (coli, mtb):
        pos = [Counter() for _ in range(NB)]
        tot_b = [0] * NB
        for _, s in org["seqs"]:
            L = len(s)
            for i, c in enumerate(s):
                b = min(NB - 1, i * NB // L)
                pos[b][c] += 1
                tot_b[b] += 1
        overall = Counter()
        overall.update(org["aa"])
        ov = sum(overall.values())
        matrix = [[round((pos[b][a] / tot_b[b]) / (overall[a] / ov), 2) if overall[a] else 0
                   for b in range(NB)] for a in AA20]
        cols = [f"{i*10}–{i*10+10}%" for i in range(NB)]
        add(f"{org['meta']['species']}: amino-acid enrichment along the chain",
            "Positional frequency divided by the proteome-wide frequency, over every real "
            "residue of every protein; 1.00 = average, darker = enriched.",
            heatmap(list(AA20), cols, matrix, color=org["meta"]["color"], unit="×", tall=1),
            wide=True)

    # ------------------------------------------------------------------ L.
    # slope chart: each of the nine reference proteomes before and after
    add("Runtime slope: InterProScan 5 to InterProScan 6",
        f"Wall time of each of the nine reference proteomes before and after; one line "
        f"per proteome, steeper = larger gain ({src}).",
        line(labels=["InterProScan 5", "InterProScan 6"],
             series=[{"name": names[i], "values": [i5[i], i6[i]]} for i in range(len(names))],
             yUnit="minutes", wide=True))

    # ------------------------------------------------------------------ M.
    # parallel coordinates: percentile profile of sampled real proteins across
    # seven bio-physical and annotation axes
    PC_METRICS = ["length", "MW", "GRAVY", "charged %", "Cys %", "matches", "IPR entries"]
    for org in (coli, mtb):
        tag = org["meta"]["species"]
        idx = {sid: i for i, (sid, _) in enumerate(org["seqs"])}
        cys = {sid: pct(s.count("C"), len(s)) for sid, s in org["seqs"]}

        def raw(sid, name):
            if name == "length":
                return org["lens"][sid]
            if name == "MW":
                return org["mw"][idx[sid]]
            if name == "GRAVY":
                return org["gravy"][idx[sid]]
            if name == "charged %":
                return org["charge"][idx[sid]] * 100
            if name == "Cys %":
                return cys[sid]
            if name == "matches":
                return org["matches"].get(sid, 0)
            return len(org["ipr_set"].get(sid, ()))

        all_ids = [sid for sid, _ in org["seqs"]]
        sv = {m: sorted(raw(sid, m) for sid in all_ids) for m in PC_METRICS}

        def pctl(sid, name):
            arr = sv[name]
            return round(100.0 * bisect_left(arr, raw(sid, name)) / (len(arr) - 1), 1)

        cand = [sid for sid in all_ids if org["matches"].get(sid)]
        sel = cand[::max(1, len(cand) // 8)][:8]
        add(f"{tag}: protein profiles across seven axes",
            "Parallel coordinates: eight sampled real proteins; each value is its "
            "percentile within the proteome. Axes, left to right: " +
            ", ".join(PC_METRICS) + ".",
            line(series=[{"name": sid, "values": [pctl(sid, m) for m in PC_METRICS]}
                         for sid in sel],
                 yMin=0, yMax=100, yUnit="percentile", wide=True))

    # ------------------------------------------------------------------ N.
    # quad-encoded bubble scatter: length x matches x InterPro entries, both proteomes
    add("Length, matches and InterPro breadth at once",
        "Each bubble is one real protein (sampled): x = length, y = matches, "
        "bubble size = distinct InterPro entries it carries.",
        scatter([{"name": org["meta"]["species"],
                  "points": sample([[org["lens"][sid], org["matches"].get(sid, 0),
                                     round(2 + 1.1 * min(8, len(org["ipr_set"].get(sid, ()))))]
                                    for sid, _ in org["seqs"] if org["matches"].get(sid)], 600),
                  "color": org["meta"]["color"]} for org in (coli, mtb)],
                xUnit="length (aa)", yUnit="matches"))

    # ------------------------------------------------------------------ O.
    # entry co-occurrence: which InterPro entries are found together in proteins
    for org in (coli, mtb):
        top12 = [e for e, _ in sorted(org["ipr_prot"].items(), key=lambda kv: -len(kv[1]))[:12]]
        sets = [org["ipr_prot"][e] for e in top12]
        mat = [[round(len(sets[i] & sets[j]) / len(sets[i] | sets[j]), 2)
                for j in range(len(top12))] for i in range(len(top12))]
        add(f"{org['meta']['species']}: InterPro entry co-occurrence",
            "Jaccard index over proteins: |A∩B| / |A∪B| for the top 12 entries by "
            "proteins; 1.00 = always seen together, 0 = never.",
            heatmap(top12, top12, mat, color=org["meta"]["color"], unit="Jaccard"),
            wide=True)

    # ------------------------------------------------------------------ P.
    ce, me = set(coli["ipr_match"]), set(mtb["ipr_match"])
    both = ce & me
    add("InterPro entries shared between the proteomes",
        f"{fmt_n(len(ce))} E. coli entries and {fmt_n(len(me))} M. tuberculosis entries; "
        "overlap = entries hit in both real annotation sets.",
        donut([{"name": "both proteomes", "value": len(both), "color": "#0b6e6e"},
               {"name": "E. coli only", "value": len(ce - me), "color": coli["meta"]["color"]},
               {"name": "M. tuberculosis only", "value": len(me - ce),
                "color": mtb["meta"]["color"]}],
              center=fmt_n(len(ce | me)), centerUnit="entries in total"))

    a0 = sum(coli["ipr_match"][e] for e in both)
    a1 = sum(coli["ipr_match"][e] for e in ce - me)
    b0 = sum(mtb["ipr_match"][e] for e in both)
    b1 = sum(mtb["ipr_match"][e] for e in me - ce)
    add("Matches from shared vs unique InterPro entries",
        "How much of each proteome's annotation traffic comes from entries the other "
        "proteome also uses.",
        stacked(["E. coli", "M. tuberculosis"],
                [{"name": "entries in both proteomes", "values": [a0, b0], "color": "#0b6e6e"},
                 {"name": "entries unique to one", "values": [a1, b1], "color": "#b4541a"}],
                yUnit="matches"))

    return charts


def main():
    charts = build()
    out = {
        "generatedBy": "scripts/build_visualizations.py",
        "sources": [
            "annotations/coli_UP000000625_annotations.tsv",
            "annotations/mtb_UP000001584_annotations.tsv",
            "data/coli_UP000000625.fasta",
            "data/mtb_UP000001584.fasta",
            "results/*.csv",
        ],
        "charts": charts,
    }
    dest = path("site", "viz-data.js")
    with open(dest, "w", encoding="utf-8") as fh:
        fh.write("// Generated by scripts/build_visualizations.py — all values come from the\n")
        fh.write("// real files listed in VIZ.sources. Do not edit by hand.\n")
        fh.write("window.VIZ=")
        json.dump(out, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write(";\n")
    size = os.path.getsize(dest)
    print(f"wrote {dest}: {len(charts)} charts, {size/1024:.0f} KB")
    assert 110 <= len(charts) <= 125, f"chart count out of range: {len(charts)}"


if __name__ == "__main__":
    main()
