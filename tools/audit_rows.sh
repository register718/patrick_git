#!/bin/bash
# Lists the reactions of Modules 3..mem whose table row has no colour mark (no source behind it).
cd "$(dirname "$0")/.."
for m in 3 4 8 9 10 11 mem; do
  awk -F' & ' -v m=$m '/&/ && !/^ID/ && $1 !~ /\\mk/ {print "M"m": "$1}' sections/tables/tab_m${m}_rxn.tex
done
