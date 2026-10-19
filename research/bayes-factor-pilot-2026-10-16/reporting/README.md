# Report production helpers

These scripts support the authored research package on `D:/heterogeneity-bayes-research-2026-10-16`. They do not replace its literature searches or regenerate reviewed prose from a model.

- `run_grok.py`: capture the explicitly requested Grok scouting/writing runs; inputs and raw output stay on D:. The writing pass disables web search. Scientific claims from either run require independent verification.
- `finalize_sources.py`: apply the first cross-review corrections to authored reports and metadata.
- `complete_package.py`: apply final editorial corrections, synchronize source ledgers, create the offline reading guide, and verify the original manuscript hashes.
- `fetch_bibliography.py`: retrieve DOI BibTeX with a cache and rate limiting. It reads final top-level reports/ledgers, excludes archived agent drafts, and records failures. Raw metadata may retain online-first years and HTML markup.
- `build_reports.py`: create self-contained HTML with MathML and PDF with Pandoc/XeLaTeX. Long prose matrices become records in the print version; original Markdown and HTML tables remain intact.
- `inspect_pdf.py`: render all PDF pages and contact sheets for manual visual inspection, plus text sanity checks. Requires the bundled Python with Pillow, pypdf, and pypdfium2.

To rerender the existing authored reports from the repository root:

```powershell
python research/bayes-factor-pilot-2026-10-16/reporting/build_reports.py --root D:/heterogeneity-bayes-research-2026-10-16 01_heterogeneity_review.md 02_bayes_factor_literature_review.md 03_standalone_research_proposal.md 04_pilot_results.md
```

Pandoc and XeLaTeX defaults point to this machine's installations and can be overridden with `--pandoc` and `--xelatex`. The root pilot README documents simulation reproduction. Preserve the raw results when changing presentation; do not replace diagnostic-selected cases in the primary operating-characteristic sample.
