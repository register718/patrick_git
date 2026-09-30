#!/usr/bin/env bash
# cited.sh sections/m3.tex   -> every key cited in that file (incl. its tables) with the PDF path.
# Bare "cited.sh" does the whole paper.  Order of the compiled numbers is irrelevant: you get keys.
cd "$(dirname "$0")/.."
files=${@:-sections/*.tex sections/tables/*.tex}
[ $# -eq 1 ] && files="$1 $(ls sections/tables/tab_$(basename ${1%.tex})_*.tex 2>/dev/null)"
cat $files | grep -o '\\cite{[^}]*}' | sed 's/\\cite{//;s/}//' | tr ',' '\n' | tr -d ' ' | sort -u | while read k; do
  f=sources/all/$k.pdf; [ -f "$f" ] && echo "$k  $f" || echo "$k  (MISSING)"; done
