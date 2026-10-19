"""Render the research Markdown as readable HTML and PDF, without changing it.

Requires pandoc and xelatex on the supplied paths. All artifacts go to --root.
The source reports and literature records are authored separately.
"""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import re
import subprocess


STYLE = """
html { color:#20252a; background:#f6f5f1; }
body { margin:0 auto; max-width:1050px; padding:4rem 3rem; background:white;
       font-family:Georgia,Cambria,serif; font-size:18px; line-height:1.65; }
h1,h2,h3,h4 { font-family:Calibri,Arial,sans-serif; line-height:1.22; color:#182d3b; }
h1 { font-size:2.35rem; margin-top:1.8em; } h2 { margin-top:2.2em; }
a { color:#166277; text-decoration-thickness:1px; text-underline-offset:3px; }
table { border-collapse:collapse; width:100%; font:15px/1.5 Calibri,Arial,sans-serif; }
th,td { border:1px solid #d4dce0; padding:10px 12px; vertical-align:top; }
th { background:#e8eef1; } tr:nth-child(even) { background:#f6f8fa; }
pre { overflow-x:auto; padding:1rem; background:#f4f6f8; }
code { font-size:.88em; } img { max-width:100%; }
#TOC { font:16px/1.6 Calibri,Arial,sans-serif; margin:2rem 0 3rem; }
blockquote { margin-left:1rem; padding-left:1rem; border-left:2px solid #b7cbd2; }
@media(max-width:700px) { body { padding:1.5rem; font-size:17px; } table { display:block; overflow-x:auto; } }
@media print { body { max-width:none; padding:0; background:white; } }
"""


def table_records_for_pdf(text: str) -> str:
    """Expand wide prose matrices into records so no unreadable tiny type is needed."""
    lines = text.splitlines()
    output, i = [], 0
    while i < len(lines):
        if (lines[i].lstrip().startswith("|") and i + 1 < len(lines)
                and re.match(r"^\s*\|[\s:|\-]+\|\s*$", lines[i + 1])):
            header = [x.strip() for x in lines[i].strip().strip("|").split("|")]
            rows, j = [], i + 2
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                rows.append([x.strip() for x in lines[j].strip().strip("|").split("|")])
                j += 1
            wide = len(header) >= 5 or any(sum(map(len, row)) > 420 for row in rows)
            if wide:
                output += ["", "*The comparison matrix is presented as individual records in this PDF; the HTML and Markdown retain the table.*", ""]
                for row in rows:
                    if len(row) != len(header):
                        raise ValueError(f"Malformed table row: {row!r}")
                    output += [f"**{row[0]}**", ""]
                    output += [f"**{h}:** {v}" + "\n" for h, v in zip(header[1:], row[1:])]
                i = j
                continue
        output.append(lines[i])
        i += 1
    return "\n".join(output) + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--pandoc", default=r"C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools\pandoc.exe")
    p.add_argument("--xelatex", default=r"C:\Users\Yusuf Balci\AppData\Roaming\TinyTeX\bin\win32\xelatex.exe")
    p.add_argument("--html-only", action="store_true")
    p.add_argument("reports", nargs="+")
    a = p.parse_args()
    build = a.root / "rendering"
    build.mkdir(parents=True, exist_ok=True)
    css = build / "report.css"
    css.write_text(STYLE, encoding="utf-8")
    header = build / "report_header.tex"
    header.write_text(r"""\usepackage{xurl}
\usepackage{microtype}
\emergencystretch=3em
\AtBeginDocument{\hypersetup{breaklinks=true}\let\setminus\backslash}
\setlength{\parskip}{5pt}
\setcounter{tocdepth}{2}
""", encoding="utf-8")
    results = []
    for name in a.reports:
        source = a.root / "reports" / name
        content = source.read_text(encoding="utf-8-sig")
        title = next((s[2:].strip() for s in content.splitlines() if s.startswith("# ")), source.stem)
        base = [a.pandoc, "-f", "markdown+tex_math_dollars+tex_math_single_backslash", "--resource-path", str(source.parent), "--standalone", "--toc", "--toc-depth=2", "-M", "lang=en-GB", "-M", "title=" + title]
        html_path = source.with_suffix(".html")
        command = base + [str(source), "--mathml", "--embed-resources", "--css", str(css), "-o", str(html_path)]
        completed = subprocess.run(command, capture_output=True, text=True)
        (build / (source.stem + "_html.log")).write_text(completed.stdout + completed.stderr, encoding="utf-8")
        completed.check_returncode()
        entry = {"source": str(source), "html": str(html_path)}
        if not a.html_only:
            pdf_input = build / (source.stem + "_print.md")
            pdf_input.write_text(table_records_for_pdf(content), encoding="utf-8")
            pdf = source.with_suffix(".pdf")
            command = base + [str(pdf_input), "--pdf-engine", a.xelatex,
                "--include-in-header", str(header), "-V", "geometry:margin=23mm",
                "-V", "papersize:a4", "-V", "fontsize:11pt", "-V", "linestretch:1.12",
                "-V", "mainfont:Cambria", "-V", "sansfont:Calibri", "-V", "monofont:Consolas",
                "-V", "colorlinks:true", "-V", "urlcolor:teal", "-V", "linkcolor:teal",
                "-o", str(pdf)]
            completed = subprocess.run(command, cwd=build, capture_output=True, text=True)
            (build / (source.stem + "_pdf.log")).write_text(completed.stdout + completed.stderr, encoding="utf-8")
            completed.check_returncode()
            entry["pdf"] = str(pdf)
        results.append(entry)
    (build / "render_manifest.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
