# Finding papers fast

Rule: **the BibTeX key is the only stable identifier.** `[17]` in the PDF changes whenever citations are added or reordered; `ford2023` does not.
Every PDF is `sources/all/<key>.pdf`; keys are `firstauthor+year` (`kae2004`, `kae2007`), so the folder is already alphabetical by author.

| I want to … | do |
|---|---|
| open a paper whose key I know / half know | `tools/open.sh ford2023` , `tools/open.sh kae` (lists if ambiguous) |
| open the paper while editing | in VS Code press `Ctrl+P`, type the key → `sources/all/<key>.pdf` |
| see all papers behind one section | `tools/cited.sh sections/m3.tex` |
| know which paper says X | `tools/find.sh "60 nM"` (full-text search, key + hit count) |
| resolve a number I see in Main.pdf | `tools/number2key.sh 17` (uses the current `Main.bbl`) or, while writing, add `\usepackage{showkeys}` temporarily |
| browse by topic | `sources/INDEX.md` (A–Z table + by-module lists) |
| use Zotero/JabRef | every `Bibliography.bib` entry has `file = {sources/all/<key>.pdf}`; import the .bib and the links work |

`sources/M*/` are convenience copies per module (same file names); edit/annotate only `sources/all/`.
