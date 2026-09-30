# Source marks: workflow

* `tools/marks/<module>.json` — one file per module (`m0` ... `m11`, `mmem`, `general`). It lists the statements of the model text ("claims")
  and, per cited paper, the passages that support them (`page`, `start` [, `end`] or an explicit `rect`).
* `python tools/assign_colors.py tools/marks/m3.json` assigns a colour (1-10, see `tools/mark_pdfs.py`) to every claim without one, such that two claims marked
  in the same paper never share a colour, and regenerates `marks_claims.tex` (claim name -> colour).
* `python tools/mark_pdfs.py tools/marks/m3.json` writes the highlighted copies to `sources_marked/<module>/<paper>.pdf` and the legend
  `sources_marked/<module>/README.md`. The originals in `sources/` are never modified. Matching ignores case, punctuation, hyphenation and understands two-column pages.
* In the model text a statement is coloured with `\mk{M3-claim-name}{text}` (running text), `\mkt{...}{text}` (table cells; no `\cite` or `\ref` inside),
  and reactions in the reaction tables with `\mkn{first}{second,third}{ID or reaction}`; the small squares are the further statements of that reaction.
  Module 0 uses the colour numbers directly (`\mk{1}{...}`), all other modules the claim names.
* `verify_m<n>.md` — result of the check of every constant, reaction and statement of the module against `dicty_reactions.jl` and the papers, including the corrections that were made.
