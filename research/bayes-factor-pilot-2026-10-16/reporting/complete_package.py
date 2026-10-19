"""Finalize editorial corrections, ledgers, and the offline reading guide."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess

from build_reports import STYLE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    root = parser.parse_args().root
    reports, sources = root / "reports", root / "sources"
    p = reports / "01_heterogeneity_review.md"
    text = p.read_text(encoding="utf-8-sig")
    paragraphs = text.split("\n\n")
    for i, paragraph in enumerate(paragraphs):
        if paragraph.startswith("[Dai and Shen (2023)]"):
            paragraphs[i] = "[Dai, Shen and Stern (2023)](https://onlinelibrary.wiley.com/doi/abs/10.1002/cjs.11728) extend the rank/U-statistic approach to observational studies using propensity-based reweighting. The final publisher citation and Crossref record list all three authors; an abbreviated author display initially obscured Stern. Reweighting changes confounding adjustment; it does not by itself turn a rank/distributional null into constant additive CATE. The article appeared online in 2022 and in the 2023 journal issue."
    p.write_text("\n\n".join(paragraphs), encoding="utf-8")
    p = reports / "02_bayes_factor_literature_review.md"
    text = p.read_text(encoding="utf-8-sig").replace("involving a prespecified modifier vector.", "involving a prespecified centered modifier vector, with \\(E[h(X)]=0\\).")
    text = text.replace(r"{\rm score}", r"{\mathrm{score}}")
    p.write_text(text, encoding="utf-8")
    p = reports / "03_standalone_research_proposal.md"
    text = p.read_text(encoding="utf-8-sig").replace("oracle is a diagnostic upper-reference calculation", "oracle is a diagnostic reference calculation")
    text = text.replace(r"{\rm eff}", r"{\mathrm{eff}}")
    marker = "**Decision after the numerical pilot.**"
    if marker not in text:
        text = text.replace("## 2. Relationship to the actual", marker + " The isolated pilot now verifies the Gaussian/subset algebra and runs it with the existing null forest. However, selected longer-chain reruns changed evidence by many log units and reversed strong null-cell evidence for heterogeneity. The current recycling estimator and short-chain fitting procedure are therefore not ready for substantive evidence claims. The next decisive experiment is a second marginal-evidence estimator for the same reduced model, together with prior and chain sensitivity. A large performance study or full horseshoe extension should wait for that numerical validation. The [completed pilot report](04_pilot_results.html) preserves the original rates and reports diagnostic-selected reruns separately; they do not establish improved calibration.\n\n## 2. Relationship to the actual", 1)
    p.write_text(text, encoding="utf-8")

    counts = {}
    for name in ["heterogeneity", "bayes"]:
        p = sources / (name + "_sources.json")
        data = json.loads(p.read_text(encoding="utf-8-sig"))
        rows = data if isinstance(data, list) else data["sources"]
        for row in rows:
            if row.get("doi") == "10.1002/cjs.11728":
                row["authors"] = "Maozhu Dai; Weining Shen; Hal S. Stern"
                row["notes"] = "Final publisher citation and Crossref agree on three authors; online 2022, issue 2023. Initial abbreviated author display was incomplete."
                row["url"] = "https://onlinelibrary.wiley.com/doi/abs/10.1002/cjs.11728"
            if row.get("doi") == "10.1016/j.spl.2024.110330":
                row["role"] = "Test-statistic Bayes factor functions and nonlocal moment priors; final publication verified and corrigendum recorded separately"
        p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with p.with_suffix(".csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        counts[name] = len(rows)
        intro = "# Final source ledger: " + ("treatment-effect heterogeneity" if name == "heterogeneity" else "Bayesian testing")
        intro += "\n\nEvidence cutoff and access date: 16 October 2026. This is a broad targeted narrative review, not a systematic-search census. Access levels distinguish primary full-text inspection, selected passages, abstracts, official documentation and inherited source audits. A linked paper is not a claim that every proof was independently audited. Blank DOI fields mean no DOI was verified in this pass.\n\nDiscovery used the current manuscript and earlier local comparator inventory, primary-paper citation chaining, subject searches, three research/code subagents and Grok CLI scouting. Searches covered global/omnibus and sparse interaction tests, orthogonal scores, kernels, qualitative effects, CATE variance, rankings, selective inference, Bayesian point-null testing, sparse priors, marginal likelihood, semiparametric theory and adjacent evidence frameworks. Grok suggestions were checked against primary sources; see `../grok/GROK_USE_AND_CORRECTIONS.md`.\n\nThe CSV/JSON and entries below incorporate both independent audits. Publication years use the final issue when established; the raw DOI cache can instead report first-online dates. Preprint and final versions are distinguished where known. Numerical experiments are documented separately and do not validate cited theorems.\n\n"
        records = []
        for i, row in enumerate(rows, 1):
            record = f"## {i}. {row['title']}\n\n{row['authors']} ({row['year']}). [Primary or official source]({row['url']})."
            if row.get("doi"):
                record += f" DOI: [{row['doi']}](https://doi.org/{row['doi']})."
            record += f"\n\n**Access:** {row['access_level']}.\n\n**Role:** {row['role']}.\n"
            if row.get("notes"):
                record += f"\n**Version or limitation:** {row['notes']}\n"
            records.append(record)
        (sources / (name + "_source_ledger.md")).write_text(intro + "\n".join(records), encoding="utf-8")

    notes = """# Bibliographic metadata and verification

`verified_references.bib` contains 60 DOI-resolved BibTeX records retrieved from Crossref or DOI content negotiation. `crossref_responses.json` preserves raw responses and `bibliography_fetch_status.csv` records retrieval status. 'Verified' here means the identifier resolved to publisher/depositor metadata, not that every bibliographic field or scientific claim was proof-audited.

These records are a metadata cache, not a fully copy-edited journal bibliography. Some depositors retain HTML tags/entities in titles or journal names, and the returned year may be the first-online year rather than the final issue year. For example, WATCH is online December 2024 / issue 2025, Zhao–Small–Ertefaie online 2021 / issue 2022, and Wakefield online 2008 / issue 2009. The narrative and source ledgers identify final years where checked. Use those records to check citation style before inserting the raw BibTeX into a manuscript.

The Dai–Shen–Stern observational paper has three authors. The final Wiley abstract page's author list and 'How to cite' section agree with Crossref; an earlier abbreviated author display omitted Stern. Datta–Guha–Shudde–Johnson's 2025 Statistics & Probability Letters paper uses DOI 10.1016/j.spl.2024.110330 and has correction DOI 10.1016/j.spl.2025.110446. Both are retained.

The two ledgers overlap, so their row counts must not be presented as the count of distinct papers. Some rows are official software documentation. Read access limitations before treating any entry as a full-text review.
"""
    (sources / "BIBLIOGRAPHY_NOTES.md").write_text(notes, encoding="utf-8")

    guide = f"""# Heterogeneity testing and Bayesian evidence

**Research package for Yusuf Balci · 16 October 2026**

The literature review and standalone proposal are complete, with an isolated code pilot using the existing stochtree forest sampler. The main practical finding is that the Gaussian evidence calculation is tractable, but its present forest-integrated estimate can be numerically unstable. The proposal now makes reliable marginal-evidence computation the next research milestone.

## Read the four reports

| Report | What it explains | Open |
|---|---|---|
| **1. Heterogeneity literature review** | The field by scientific null, main competitors, assumptions, blind spots, a comparison matrix and your paper's positioning | [PDF](reports/01_heterogeneity_review.pdf) · [HTML](reports/01_heterogeneity_review.html) · [Markdown](reports/01_heterogeneity_review.md) |
| **2. Bayes-factor literature review** | Held–Ott, point versus interval nulls, sparse priors, Savage–Dickey, statistic-based evidence, semiparametric theory and numerical methods | [PDF](reports/02_bayes_factor_literature_review.pdf) · [HTML](reports/02_bayes_factor_literature_review.html) · [Markdown](reports/02_bayes_factor_literature_review.md) |
| **3. Standalone research proposal** | A precise model, evidence derivations, candidate theorem, fair experiments, work packages and stopping criteria; extensions to ATE and other restrictions | [PDF](reports/03_standalone_research_proposal.pdf) · [HTML](reports/03_standalone_research_proposal.html) · [Markdown](reports/03_standalone_research_proposal.md) |
| **4. Completed pilot results** | Actual simulations, uncertainty intervals, numerical failures, prior sensitivity and longer-chain stress checks | [PDF](reports/04_pilot_results.pdf) · [HTML](reports/04_pilot_results.html) · [Markdown](reports/04_pilot_results.md) |

For a first pass, read this page, sections 1 and 12–14 of the field review, sections 2–3 and 7 of the Bayesian review, and the pilot's summary. Then use the proposal as the agenda for a supervisor meeting.

## The scientific picture

**Heterogeneity is several different questions.** A patient's baseline prognosis can vary while the treatment's added benefit stays constant. A CATE homogeneity test asks about that added mean benefit. A finite-basis test asks whether selected directions show variation; a sign-change test asks whether treatment helps some and harms others; a ranking or policy test asks whether personalization is useful. These are related targets, not interchangeable competitors. The field review explains which comparisons are fair.

**The closest prior art is broader than BCF.** The main testing comparisons include fixed-basis orthogonal projections, Qian et al.'s sparse semiparametric interaction tests, kernel score methods, and recent omnibus/rule-process tests. Semi-BART is close model-architecture prior art. GenericML, GRF calibration, RATE, CATE-variance methods and selective inference address important neighboring questions. Merely combining a forest nuisance with linear modifiers is not a sufficient novelty claim. See the cited papers and qualification of each result in [report 1](reports/01_heterogeneity_review.html).

**A Bayes factor compares two specified models.** Posterior model odds equal prior odds multiplied by the Bayes factor. Evidence favoring a constant-effect model is relative to the alternative prior actually used; it is not proof of universal homogeneity. A continuous horseshoe prior alone assigns no positive probability to an exact coefficient-zero model. A global model indicator can supply that point mass. The pilot instead starts with proper Gaussian slabs and a sparse nonempty-subset mixture, which can be integrated analytically. [Report 2](reports/02_bayes_factor_literature_review.html) explains the supervisor's Held–Ott reference and the surrounding literature.

**The first BF pilot exposes a computational problem worth solving.** The algebra passed direct-density and independent quadrature checks, but some short-chain forest estimates were dominated by only a few draws. In a deliberately selected confounded null dataset, dense log BF changed from about 20.17 to −0.47 in the longer run; a randomized null case changed from 12.54 to −2.78. Those are changes in the natural logarithm of the BF and reverse strong apparent evidence for heterogeneity. These reruns are diagnostics, not an independent performance sample or proof that the longer estimates are correct. [Report 4](reports/04_pilot_results.html) retains the original simulation results and reports these cases separately.

**The existing score remains an essential benchmark.** A fixed-dimension Gaussian-score BF is monotone in the current quadratic score for a fixed prior parameter. It cannot improve ranking-based detection at matched thresholds. Full-model sparse evidence can behave differently, but unequal conventional cutoffs, unequal prior signal energy and unstable integration prevent a superiority claim from the present pilot.

## Recommended next research step

Keep the Gaussian, fixed-basis, known-propensity model as the reference experiment and obtain an independent marginal-evidence calculation, such as a carefully derived bridge or product-space method. Compare the two calculations across null and alternative cases and across scientific prior scales. Then run a separately calibrated, larger comparison with the heterogeneity score. Extend to the actual linked/unlinked horseshoe hierarchy only after the reference calculation is reliable. The proposal sets out what needs proving and what would make the project worth continuing.

The pilot uses fixed known noise variance and fixed leaf variance; it does not establish a Bayes factor for every option in `bcf_linear_probit`, binary outcomes, estimated propensity, or the manuscript's exact practical defaults. Both ordinary modifiers and modifier–modifier products must be inactive for a constant-CATE null. Treatment centering and compatible nuisance priors also form part of the model definition.

## Sources, code and reproducibility

There are **{counts['heterogeneity']} annotated field-review records** and **{counts['bayes']} Bayesian-review records**, with overlap, plus **60 DOI-resolved BibTeX records**. This is a broad targeted review, not a registered systematic review or a claim to cover every paper. Source access and preprint status are recorded explicitly.

- [Heterogeneity source ledger](sources/heterogeneity_source_ledger.md) · [CSV](sources/heterogeneity_sources.csv)
- [Bayesian source ledger](sources/bayes_source_ledger.md) · [CSV](sources/bayes_sources.csv)
- [Retrieved BibTeX](sources/verified_references.bib) · [metadata limitations](sources/BIBLIOGRAPHY_NOTES.md)
- [Independent review of the field map](sources/REVIEW_HETEROGENEITY.md) · [independent review of Bayesian derivations and pilot](sources/REVIEW_BAYES_AND_PILOT.md)
- [Grok research/writing use and corrections](grok/GROK_USE_AND_CORRECTIONS.md). Raw prompts and responses are preserved; unchecked suggestions are not scientific references.
- [Pilot results folder](pilot-results/) contains per-replication data, summaries, diagnostics, plots, seeds, settings, runtime and provenance.

Code remains at `C:/Users/Yusuf Balci/Documents/stochtree/research/bayes-factor-pilot-2026-10-16/`. Its `README.md` gives executable reproduction commands and the exact scope. No existing package source, exports, native implementation, manuscript, or existing tests were changed. New functionality is opt-in research code; no commit or release was made.

The PDFs are designed for reading and sharing; HTML files retain wide comparison matrices and can be opened offline. Markdown remains editable. All research outputs and large simulation artifacts are on D:.
"""
    (root / "START_HERE.md").write_text(guide, encoding="utf-8")
    pandoc = r"C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools\pandoc.exe"
    css = root / "rendering/report.css"
    css.write_text(STYLE, encoding="utf-8")
    subprocess.run([pandoc, str(root / "START_HERE.md"), "--standalone", "--embed-resources", "--css", str(css), "-M", "title=Heterogeneity testing and Bayesian evidence", "-o", str(root / "index.html")], check=True)
    before = json.loads((sources / "input_manuscript_hashes.json").read_text(encoding="utf-8-sig"))
    checked = []
    for row in before:
        now = hashlib.sha256(Path(row["Path"]).read_bytes()).hexdigest().upper()
        checked.append(dict(path=row["Path"], original_sha256=row["Hash"], final_sha256=now, unchanged=now == row["Hash"]))
    (sources / "manuscript_unchanged_check.json").write_text(json.dumps(checked, indent=2), encoding="utf-8")
    assert all(r["unchanged"] for r in checked), "Manuscript changed during research"
    print(json.dumps({"source_records": counts, "manuscripts_unchanged": len(checked), "index": str(root / "index.html")}, indent=2))


if __name__ == "__main__":
    main()
