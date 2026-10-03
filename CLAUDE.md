# Vorgehen: Aussage ohne Quelle belegen (cite + Farbmarkierung + Markierung im PDF)

Details der Tools: `tools/README.md`, `tools/marks/README.md`.

1. Aussage in `sections/*.tex` (oder `sections/tables/`) finden; prüfen, welche Papers in `sources/all/<key>.pdf` (BibTeX-Key = Dateiname) sie belegen
   (`tools/find.sh "text"`, oder `pdftotext`/pymupdf-Volltextsuche über `sources/all`). Nur Passagen zitieren, die die Aussage wörtlich stützen; wenn kein Paper vorhanden ist, `sources/TO_DOWNLOAD.md` ergänzen statt zu raten.
2. In `tools/marks/<modul>.json` einen Claim (`claims.<id>` mit `text`, `color`) anlegen und je Paper unter `papers.<key>` eine Markierung (`claim`, `page`, `start`, optional `end`) ergänzen. Farbe 1-10 wählen, die in *diesem Paper* noch nicht für einen anderen Claim benutzt ist (Paper-Liste der Marks prüfen).
3. Im TeX: `\mk{<Modul>-<claim-id mit Bindestrichen>}{Aussage}~\cite{key1,key2}` (Zitat direkt hinter der markierten Aussage). In Reaktionstabellen weitere Farben als Quadrate: `\mkn{erste}{zweite,dritte}{ID}`.
4. `python3 tools/assign_colors.py tools/marks/<modul>.json` (schreibt `marks_claims.tex`), dann `python3 tools/mark_pdfs.py tools/marks/<modul>.json` (schreibt `sources_marked/<modul>/*.pdf` + README); in der Ausgabe darf kein "not found" stehen.
5. `latexmk -pdf Main.tex`, auf undefined citations/references prüfen.

## Prüfung: hat jede Reaktion (ab Modul 3, Reaktion 3.3) eine Quelle?

`tools/audit_rows.sh` listet alle Zeilen der Reaktionstabellen (M3, M4, M8-M11, Mmem) ohne Farbmarkierung. Stand 2026-10-03: nur 11.D4a-d und 11.E2 (Abbau/Turnover von PdsA/PdiA) sind ohne Quelle; dafür gibt es keine Messung, das ist im Text als Designentscheidung gekennzeichnet und steht in `sources/TO_DOWNLOAD.md`. Auch 8.9/8b.3 (PKB phosphoryliert CRAC) ist im Text als Annahme gekennzeichnet; belegt ist nur, dass TORC2, PKBs und CRAC für ACA nötig sind.
Ehrlichkeit vor Vollständigkeit: wenn kein Paper die Aussage trägt, nicht eine unpassende Stelle markieren, sondern als Annahme/Designentscheidung kennzeichnen. Einen Passagen-Text immer im PDF prüfen (OCR-Fehler, z. B. franke1981: "10 µM to 2 mM").
