# Insert the generated Husler-Reiss table bodies into the manuscript.

import re
import pathlib

root = pathlib.Path(__file__).resolve().parent.parent
tex = root / "pgm.tex"
gen = (root / "results" / "hr_tables.tex").read_text(encoding="utf-8")

blocks, cur = {}, None
for line in gen.splitlines():
    m = re.match(r"%% ---- TABLE (\S+)", line)
    if m:
        cur = m.group(1).rstrip(")").rstrip("(")
        cur = m.group(1)
        blocks[cur] = []
        continue
    if cur and line.strip():
        blocks[cur].append(line)

t2 = blocks["2"]
t3 = blocks["3"]
t4a = blocks["4(a)"]
t4b = blocks["4(b)"]
t4c = blocks["4(c)"]

src = tex.read_text(encoding="utf-8")

def replace_between(text, start_marker, end_marker, new_body, label):
    """Replace everything strictly between start_marker and end_marker."""
    i = text.index(start_marker) + len(start_marker)
    j = text.index(end_marker, i)
    print(f"  {label}: replaced {text[i:j].count(chr(92)+chr(92))} rows")
    return text[:i] + new_body + text[j:]

def cells(row):
    return [c.strip() for c in row.rstrip("\\ ").split("&")]

e, nb, na, gl, sb, sa = (cells(t2[0]), cells(t2[1]), cells(t2[2]),
                         cells(t2[3]), cells(t2[5]), cells(t2[6]))
table2 = f"""
\\EMTPtwo{{}}                         & (none) & --        & {e[3]}  & {e[4]}           & {e[5]} & {e[6]} \\\\
\\midrule
\\multirow{{2}}{{*}}{{eglearn-NS}}        & BIC    & {nb[2]}   & {nb[3]}  & {nb[4]}  & {nb[5]} & \\multirow{{2}}{{*}}{{{nb[6]}}} \\\\
                                   & AIC    & {na[2]}   & {na[3]} & {na[4]}           & {na[5]} & \\\\
\\midrule
eglearn-glasso                     & BIC \\& AIC & {gl[2]} & {gl[3]} & {gl[4]}         & {gl[5]} & {gl[6]} \\\\
\\midrule
\\multirow{{2}}{{*}}{{Spectral-HR}}       & BIC    & {sb[2]}       & {sb[3]} & {sb[4]}   & {sb[5]} & \\multirow{{2}}{{*}}{{{sb[6]}}} \\\\
                                   & AIC    & {sa[2]}      & {sa[3]}  & {sa[4]}           & {sa[5]} & \\\\
"""
src = replace_between(
    src,
    "Method & Selected by & Param & Edges & BIC & AIC & Time (s) \\\\\n\\midrule\n",
    "\\bottomrule",
    table2.lstrip("\n"),
    "Table 2",
)

table3 = "\n".join(t3) + "\n"
src = replace_between(
    src,
    "Method & Selected by & Param & Edges & Train BIC & Train AIC "
    "& HR test log-lik.\\ & Time (s) \\\\\n\\midrule\n",
    "\\bottomrule",
    table3,
    "Table 3",
)

pos = 0
for name, body in (("4(a)", t4a), ("4(b)", t4b), ("4(c)", t4c)):
    header = "$\\eta$ & Edges & BIC & AIC \\\\\n\\midrule\n" if name == "4(a)" \
        else "$\\rho$ & Edges & BIC & AIC \\\\\n\\midrule\n"
    i = src.index(header, pos) + len(header)
    j = src.index("\\bottomrule", i)
    new_body = "\n".join(body) + "\n"
    src = src[:i] + new_body + src[j:]
    pos = i + len(new_body)
    print(f"  Table {name}: {len(body)} rows")

tex.write_text(src, encoding="utf-8")
print("pgm.tex updated")
