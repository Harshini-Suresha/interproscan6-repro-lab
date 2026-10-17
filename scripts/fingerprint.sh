#!/usr/bin/env bash
# fingerprint.sh output.tsv  -> hash of sorted annotations with the date column removed.
# Columns assume standard InterProScan TSV (date is column 11). Check your file first.
cut -f1-10,12-15 "$1" | sort | sha256sum | cut -c1-12
