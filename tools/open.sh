#!/usr/bin/env bash
# open.sh <bibtex-key | number-free fragment>   e.g.  ./tools/open.sh ford2023   or   ./tools/open.sh kae
# Opens sources/all/<key>.pdf ; with a fragment that matches several keys, lists them.
cd "$(dirname "$0")/.."
m=( $(ls sources/all/ | grep -i -- "$1" | grep -v _supp) )
if   [ ${#m[@]} -eq 1 ]; then xdg-open "sources/all/${m[0]}" >/dev/null 2>&1 &
elif [ ${#m[@]} -eq 0 ]; then echo "no paper matches '$1'"
else printf '%s\n' "${m[@]%.pdf}"; fi
