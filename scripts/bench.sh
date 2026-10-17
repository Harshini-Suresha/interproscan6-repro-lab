#!/usr/bin/env bash
# bench.sh NAME [extra nextflow args]
# Follows the official README for InterProScan 6.0.1. Put your FASTA at data/proteins.faa first.
name=$1; shift; mkdir -p runs/$name; cd runs/$name
/usr/bin/time -v -o time.txt \
  nextflow run ebi-pf-team/interproscan6 -r 6.0.1 -profile docker \
  --input ../../data/proteins.faa --datadir ../../data --interpro 108.0 "$@"
grep -E "Elapsed|Maximum resident" time.txt
