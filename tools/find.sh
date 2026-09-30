#!/usr/bin/env bash
# find.sh "phrase"  -> which papers contain the phrase (full text, case-insensitive), by BibTeX key.
# Use it to find "where did I read that K_d is 60 nM".  Text is cached in tools/.txt after the first run.
cd "$(dirname "$0")/.."; mkdir -p tools/.txt
for f in sources/all/*.pdf; do t=tools/.txt/$(basename "${f%.pdf}").txt
  [ "$t" -nt "$f" ] || pdftotext -q "$f" "$t" 2>/dev/null
done
grep -i -c -- "$1" tools/.txt/*.txt | grep -v ':0$' | sed 's|tools/.txt/||;s|\.txt:| : |' | sort -t: -k2 -nr
