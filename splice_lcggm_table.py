# Insert the generated LCGGM table body into the manuscript.

import pathlib

root = pathlib.Path(__file__).resolve().parent.parent
tex = root / "pgm.tex"
gen = (root / "results" / "lcggm_table.tex").read_text(encoding="utf-8")

body = "\n".join(l for l in gen.splitlines() if not l.startswith("%%")).strip("\n") + "\n"

HEADER = ("Setting & Method & Edges & Precision & Recall & F1 & Test log-lik.\\ "
          "& Time (s) \\\\\n\\midrule\n")

src = tex.read_text(encoding="utf-8")
if HEADER not in src:
    raise SystemExit("Table 1 header not found in pgm.tex -- did the columns change?")
i = src.index(HEADER) + len(HEADER)
j = src.index("\\bottomrule", i)
old_rows = src[i:j].count("\\\\")
src = src[:i] + body + src[j:]
tex.write_text(src, encoding="utf-8", newline="\n")
print(f"Table 1: replaced {old_rows} rows with {body.count(chr(92) + chr(92))}")
