#!/usr/bin/env bash
# number2key.sh 17   -> the key behind [17] in the CURRENT compiled Main.pdf (reads Main.bbl, so recompile first).
cd "$(dirname "$0")/.."
awk -v n="$1" '/\\bibitem/ {c++} c==n && /\\bibitem/ {print}' Main.bbl | head -1
