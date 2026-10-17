#!/usr/bin/env python3
"""fetch_real_data.py

Download two real UniProt reference proteomes (the paper's two smallest) and
their real InterPro member-database cross-references, then write InterProScan
style annotation files and summary CSVs.

Sources
  Sequences : UniProt REST, stream endpoint (proteome:UP000000625 / UP000001584)
  Annotations: UniProt REST search, the xref_* TSV fields (Pfam, PROSITE, SMART,
               SUPERFAMILY, PRINTS, PIRSF, CATH-Gene3D, CDD, HAMAP, NCBIfam,
               PANTHER, SFLD, InterPro)
  Names / signature -> InterPro mapping: InterPro interpro.xml.gz (release file)

Nothing here runs InterProScan. The annotations are the cross-references that
UniProt already publishes, reformatted to InterProScan's 15-column TSV so the
site's Compare page and scripts/compare.py work unchanged.

Usage: python3 scripts/fetch_real_data.py
Writes: data/, annotations/, results/reference_summary.csv,
        results/matches_by_signature.csv, results/real_data_sources.csv
"""
import gzip
import json
import os
import re
import sys
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
ANN = os.path.join(ROOT, "annotations")
RES = os.path.join(ROOT, "results")
CACHE = os.path.join(DATA, ".cache")

UNIPROT = "https://rest.uniprot.org/uniprotkb"
INTERPRO_XML = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz"

PROTEOMES = [
    ("UP000000625", "coli", "Escherichia coli", "E. coli"),
    ("UP000001584", "mtb", "Mycobacterium tuberculosis", "M. tuberculosis"),
]

# UniProt xref field -> (analysis label, db attr in interpro.xml, id filter)
XREFS = [
    ("xref_interpro", "InterPro", None, None),
    ("xref_pfam", "Pfam", "PFAM", None),
    ("xref_prints", "PRINTS", "PRINTS", None),
    ("xref_prosite", None, None, "prosite"),
    ("xref_smart", "SMART", "SMART", None),
    ("xref_supfam", "SUPERFAMILY", "SSF", None),
    ("xref_pirsf", "PIRSF", "PIRSF", None),
    ("xref_gene3d", "CATH-Gene3D", "CATHGENE3D", None),
    ("xref_cdd", "CDD", "CDD", None),
    ("xref_hamap", "HAMAP", "HAMAP", None),
    ("xref_ncbifam", "NCBIfam", "NCBIFAM", None),
    ("xref_panther", "PANTHER", "PANTHER", None),
    ("xref_sfld", "SFLD", "SFLD", None),
]
FIELDS = ["accession", "length"] + [x[0] for x in XREFS]


def http_get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "interproscan6-repro-lab/1.0"})
            return urllib.request.urlopen(req, timeout=120)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and i < tries - 1:
                time.sleep(3 * (i + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if i < tries - 1:
                time.sleep(3 * (i + 1))
                continue
            raise


def next_url(headers):
    link = headers.get("Link", "")
    m = re.search(r'<([^>]+)>;\s*rel="next"', link)
    return m.group(1) if m else None


def download_fasta(proteome, out):
    url = f"{UNIPROT}/stream?query=proteome:{proteome}&format=fasta"
    r = http_get(url)
    release = r.headers.get("X-UniProt-Release", "")
    release_date = r.headers.get("X-UniProt-Release-Date", "")
    with open(out, "wb") as fh:
        fh.write(r.read())
    n = sum(1 for line in open(out) if line.startswith(">"))
    return n, release, release_date


def fetch_tsv(proteome):
    """Return (rows, header, release, release_date); rows are lists of cells."""
    url = (f"{UNIPROT}/search?query=proteome:{proteome}&format=tsv&size=500"
           f"&fields={','.join(FIELDS)}")
    rows, header = [], None
    release = release_date = ""
    while url:
        r = http_get(url)
        release = r.headers.get("X-UniProt-Release", release) or release
        release_date = r.headers.get("X-UniProt-Release-Date", release_date) or release_date
        text = r.read().decode("utf-8")
        lines = text.splitlines()
        if header is None:
            header = lines[0].split("\t")
            lines = lines[1:]
        for line in lines:
            if line.strip():
                rows.append(line.split("\t"))
        url = next_url(r.headers)
    return rows, header, release, release_date


def ensure_interpro_xml():
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "interpro.xml.gz")
    if not os.path.exists(path):
        sys.stderr.write("downloading interpro.xml.gz (~43 MB)...\n")
        r = http_get(INTERPRO_XML)
        with open(path, "wb") as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    return path


def load_interpro_names(path):
    """Stream interpro.xml.gz -> (ipr_name, sig_to_ipr, sig_name, release)."""
    ipr_name, sig_to_ipr, sig_name = {}, {}, {}
    release = ""
    cur_ipr = cur_name = None
    in_member = False
    ipr_text_parts = []
    with gzip.open(path, "rb") as fh:
        for event, elem in ET.iterparse(fh, events=("start", "end")):
            tag = elem.tag
            if event == "start":
                if tag == "interpro":
                    cur_ipr = elem.get("id")
                    ipr_text_parts = []
                elif tag == "member_list":
                    in_member = True
                elif tag == "db_xref" and in_member and cur_ipr:
                    dbkey = elem.get("dbkey")
                    name = elem.get("name") or ""
                    if dbkey:
                        sig_to_ipr.setdefault(dbkey, cur_ipr)
                        if name:
                            sig_name.setdefault(dbkey, name)
                elif tag == "dbinfo" and elem.get("dbname") == "INTERPRO":
                    release = elem.get("version", "")
            else:
                if tag == "name" and cur_ipr and not in_member and not ipr_text_parts:
                    ipr_text_parts.append((elem.text or "").strip())
                elif tag == "member_list":
                    in_member = False
                elif tag == "interpro":
                    if cur_ipr and ipr_text_parts:
                        ipr_name[cur_ipr] = ipr_text_parts[0].replace("\n", " ").strip()
                    cur_ipr = None
                    elem.clear()
    return ipr_name, sig_to_ipr, sig_name, release


def prosite_analysis(acc):
    return "PROSITE profiles" if acc.startswith("PS5") else "PROSITE Patterns"


def bare_accession(header_token):
    """UniProt FASTA tokens are 'sp|ACC|NAME'; return ACC, else the token."""
    parts = header_token.split("|")
    return parts[1] if len(parts) == 3 else parts[0]


def fasta_accessions(path):
    accs = set()
    for line in open(path):
        if line.startswith(">"):
            accs.add(bare_accession(line[1:].strip().split()[0]))
    return accs


def build_annotations(rows, header, ipr_name, sig_to_ipr, sig_name, keep):
    """Return (records, annotated) where records are TSV rows.

    UniProt returns TSV columns in the order requested, so fields are matched
    positionally against FIELDS rather than by the human-readable header.
    Only proteins in `keep` (the FASTA accession set) are written, because the
    search endpoint can return a few entries beyond the streamed proteome.
    """
    if len(header) != len(FIELDS):
        raise SystemExit(f"unexpected TSV columns: {header}")
    idx = {field: i for i, field in enumerate(FIELDS)}
    records = []
    annotated = set()
    date = "-"
    for row in rows:
        acc = row[idx["accession"]]
        if acc not in keep:
            continue
        length = row[idx["length"]]
        for field, label, dbattr, special in XREFS:
            cell = row[idx[field]] if field in idx else ""
            if not cell or cell == "-":
                continue
            for sig in [s for s in cell.split(";") if s]:
                if label == "InterPro":
                    analysis = "InterPro"
                    desc = ipr_name.get(sig, "-")
                    ipr = sig
                    ipr_desc = desc
                else:
                    if special == "prosite":
                        analysis = prosite_analysis(sig)
                    else:
                        analysis = label
                    desc = sig_name.get(sig, "-")
                    ipr = sig_to_ipr.get(sig, "-")
                    ipr_desc = ipr_name.get(ipr, "-") if ipr != "-" else "-"
                records.append("\t".join([
                    acc, "-", length, analysis, sig, desc,
                    "-", "-", "-", "T", date, ipr, ipr_desc, "-", "-",
                ]))
                annotated.add(acc)
    records.sort()
    return records, annotated


def main():
    for d in (DATA, ANN, RES):
        os.makedirs(d, exist_ok=True)
    xml_path = ensure_interpro_xml()
    sys.stderr.write("parsing interpro.xml.gz...\n")
    ipr_name, sig_to_ipr, sig_name, ipr_release = load_interpro_names(xml_path)
    sys.stderr.write(f"InterPro {ipr_release}: {len(ipr_name)} entries, "
                     f"{len(sig_to_ipr)} signatures\n")

    sources = []
    summary = []
    matches = []
    for proteome, slug, species, short in PROTEOMES:
        fasta = os.path.join(DATA, f"{slug}_{proteome}.fasta")
        sys.stderr.write(f"downloading {short} FASTA...\n")
        n_fasta, release, release_date = download_fasta(proteome, fasta)
        sys.stderr.write(f"  {n_fasta} sequences (UniProt {release})\n")
        keep = fasta_accessions(fasta)
        sys.stderr.write(f"fetching {short} annotations...\n")
        rows, header, release, release_date = fetch_tsv(proteome)
        records, annotated = build_annotations(rows, header, ipr_name, sig_to_ipr, sig_name, keep)
        annpath = os.path.join(ANN, f"{slug}_{proteome}_annotations.tsv")
        with open(annpath, "w") as fh:
            fh.write("\n".join(records) + "\n")

        # per-dataset summaries (denominator = streamed proteome size)
        n = n_fasta
        sigs = {r.split("\t")[4] for r in records if r.split("\t")[3] != "InterPro"}
        iprs = {r.split("\t")[11] for r in records if r.split("\t")[11] != "-"}
        interpro_annotated = {
            r.split("\t")[0] for r in records if r.split("\t")[3] == "InterPro"
        }
        summary.append({
            "dataset": slug, "species": species, "uniprot_proteome": proteome,
            "n_proteins": n, "proteins_with_match": len(annotated),
            "coverage_pct": round(100.0 * len(annotated) / n, 2) if n else 0,
            "proteins_with_interpro": len(interpro_annotated),
            "interpro_coverage_pct": round(100.0 * len(interpro_annotated) / n, 2) if n else 0,
            "total_matches": len(records), "distinct_signatures": len(sigs),
            "distinct_interpro_entries": len(iprs),
            "annotation_file": os.path.relpath(annpath, ROOT),
        })
        sources.append({
            "dataset": slug, "species": species, "uniprot_proteome": proteome,
            "sequences": n_fasta, "uniprot_release": release,
            "uniprot_release_date": release_date, "interpro_release": ipr_release,
            "fasta": os.path.relpath(fasta, ROOT), "source": "UniProt REST + InterPro",
        })
        # matches by signature (member databases, not the InterPro rollup)
        counts = {}
        for r in records:
            f = r.split("\t")
            if f[3] == "InterPro":
                continue
            key = (f[3], f[4], f[5], f[11], f[12])
            c = counts.setdefault(key, [0, set()])
            c[0] += 1
            c[1].add(f[0])
        for (db, sig, desc, ipr, iprd), (m, prots) in counts.items():
            matches.append({
                "dataset": slug, "member_database": db, "signature": sig,
                "description": desc, "interpro": ipr,
                "interpro_description": iprd, "matches": m,
                "proteins_with_hit": len(prots),
            })
        sys.stderr.write(f"  wrote {len(records)} matches for {slug}\n")

    matches.sort(key=lambda x: (x["dataset"], -x["matches"]))
    _write_csv(os.path.join(RES, "reference_summary.csv"), summary)
    _write_csv(os.path.join(RES, "matches_by_signature.csv"), matches)
    _write_csv(os.path.join(RES, "real_data_sources.csv"), sources)
    sys.stderr.write("done\n")


def _write_csv(path, rows):
    if not rows:
        return
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
