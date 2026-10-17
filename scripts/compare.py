#!/usr/bin/env python3
"""compare.py A.tsv B.tsv [n_proteins] [sig|ipr]  -> coverage and agreement (Jaccard)."""
import sys
def load(path, col):
    pairs, prots = set(), set()
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or line.startswith("#"): continue
        prots.add(f[0])
        if f[col] not in ("", "-", "NULL"): pairs.add((f[0], f[col]))
    return pairs
a_path, b_path = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else None
col = 11 if (len(sys.argv) > 4 and sys.argv[4] == "ipr") else 4
A, B = load(a_path, col), load(b_path, col)
cov = lambda s: len({p for p, _ in s})
n = n or max(cov(A), cov(B))
u = len(A | B)
print(f"coverage A {100*cov(A)/n:.1f}%  coverage B {100*cov(B)/n:.1f}%")
print(f"shared {len(A&B)}  only A {len(A-B)}  only B {len(B-A)}")
print(f"agreement (Jaccard) {100*len(A&B)/u:.1f}%" if u else "no annotations")
