#!/usr/bin/env python3
"""Highlight the passages of cited papers that support a statement of the model text.

Usage:  python tools/mark_pdfs.py tools/marks/m0.json

Spec format (see tools/marks/README.md): every mark names a paper, a PDF page (1-based), the
first words of the passage (`start`) and optionally its last words (`end`); everything from
the first word of `start` to the last word of `end` is highlighted in the colour of the claim.
Matching ignores case, punctuation, line breaks and hyphenation.

Output: sources_marked/<module>/<paper>.pdf (the originals in sources/ are never touched).
"""
import json
import re
import sys
import unicodedata
from pathlib import Path

import pymupdf

ROOT = Path(__file__).resolve().parent.parent

# colour id -> (RGB 0..255).  The same table is used for the LaTeX macros in Main.tex.
PALETTE = {
    1: (255, 235, 59),   # yellow
    2: (139, 220, 110),  # green
    3: (100, 215, 240),  # cyan
    4: (255, 150, 200),  # pink
    5: (255, 178, 80),   # orange
    6: (200, 160, 255),  # violet
    7: (215, 180, 130),  # tan
    8: (140, 170, 255),  # blue
    9: (255, 120, 110),  # salmon
    10: (190, 190, 190), # grey
}


def norm(s):
    s = unicodedata.normalize("NFKD", s)
    return re.sub(r"[^0-9a-z]", "", s.lower())


def _split(page, words):
    """x position of the column gap (fewest words crossing it) in the middle 20 % of the page, or None."""
    w = page.rect.width
    best = None
    for x in range(int(0.4 * w), int(0.6 * w), 2):
        crossing = sum(1 for t in words if t[0] < x < t[2])
        if best is None or crossing < best[0]:
            best = (crossing, x)
    return best[1] if best and best[0] <= 0.02 * len(words) else None


def _orders(page):
    """Candidate reading orders of the words of a page: as stored, and two-column (left column first)."""
    words = page.get_text("words")
    yield words
    x = _split(page, words)
    key = lambda t: (round(t[1] / 3), t[0])
    if x:
        left = [t for t in words if (t[0] + t[2]) / 2 < x]
        right = [t for t in words if (t[0] + t[2]) / 2 >= x]
        yield sorted(left, key=key) + sorted(right, key=key)
    yield sorted(words, key=key)


def locate(page, start, end):
    """Return the list of word tuples between the first match of start and the end of end."""
    s, e = norm(start), norm(end) if end else None
    last = None
    for words in _orders(page):
        joined, spans = "", []
        for w in words:
            n = norm(w[4])
            spans.append((len(joined), len(joined) + len(n)))
            joined += n
        i = joined.find(s)
        if i < 0:
            last = f"start not found on page {page.number + 1}: {start!r}"
            continue
        j = i + len(s)
        if e:
            k = joined.find(e, i)
            if k < 0:
                last = f"end not found on page {page.number + 1}: {end!r}"
                continue
            j = k + len(e)
        if joined.find(s, i + 1) >= 0:
            print(f"    warning: start is ambiguous on page {page.number + 1}: {start[:50]!r}")
        return [w for w, (a, b) in zip(words, spans) if b > i and a < j and b > a]
    raise ValueError(last)


def line_rects(words):
    """One rectangle per text line (merged words), so that the highlight is contiguous."""
    lines = {}
    for w in words:
        key = (w[5], w[6])
        r = pymupdf.Rect(w[:4])
        lines[key] = lines[key] | r if key in lines else r
    return list(lines.values())


def mark(spec_path):
    spec = json.loads(Path(spec_path).read_text())
    module = spec["module"]
    src_dir = ROOT / "sources" / module
    out_dir = ROOT / "sources_marked" / module
    out_dir.mkdir(parents=True, exist_ok=True)
    claims = spec["claims"]
    legend = []
    for paper, marks in spec["papers"].items():
        doc = pymupdf.open(src_dir / f"{paper}.pdf")
        print(f"{paper}: {len(marks)} marks")
        for m in marks:
            cid = m["claim"]
            color = PALETTE[claims[cid]["color"]]
            page = doc[m["page"] - 1]
            if "rect" in m:  # equations and figures: explicit rectangle in PDF points
                rects, words = [pymupdf.Rect(m["rect"])], [(0, 0, 0, 0, m.get("label", "<rect>"))]
            else:
                words = locate(page, m["start"], m.get("end"))
                rects = line_rects(words)
            rgb = [c / 255 for c in color]
            if "rect" in m:
                annot = page.add_rect_annot(rects[0])
                annot.set_colors(stroke=rgb, fill=rgb)
                annot.set_border(width=0.5)
                annot.set_opacity(0.35)
            else:
                annot = page.add_highlight_annot(quads=[r.quad for r in rects])
                annot.set_colors(stroke=rgb)
                annot.set_opacity(0.55)
            annot.set_info(title=f"{module} {cid}", content=f"{cid}: {claims[cid]['text']}")
            annot.update()
            legend.append((paper, m["page"], claims[cid]["color"], cid, claims[cid]["text"]))
            text = " ".join(w[4] for w in words)
            print(f"  [{cid}] p{m['page']}: {text[:110]}{'...' if len(text) > 110 else ''}")
        doc.save(out_dir / f"{paper}.pdf", garbage=3, deflate=True)
        doc.close()
    write_legend(out_dir / "README.md", module, legend)


def write_legend(path, module, legend):
    names = {1: "yellow", 2: "green", 3: "cyan", 4: "pink", 5: "orange", 6: "violet",
             7: "tan", 8: "blue", 9: "salmon", 10: "grey"}
    lines = [f"# {module}: marked passages", "",
             f"Colour = statement of `sections/*.tex` (highlighted in the same colour via `\\mk`/`\\mkt`). "
             "Generated by `tools/mark_pdfs.py`; do not edit by hand.", "",
             "| Paper | Page | Colour | Id | Statement supported by the marked passage |", "|---|---|---|---|---|"]
    for paper, page, col, cid, text in sorted(legend, key=lambda x: (x[0], x[1])):
        lines.append(f"| {paper} | {page} | {names[col]} | {cid} | {text} |")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        mark(p)
