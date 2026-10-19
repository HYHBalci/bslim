"""Apply independently checked literature-audit corrections to the D: reports.

Idempotent editorial updates. Does not modify manuscript or package code.
"""
import argparse
import csv
import json
from pathlib import Path


def insert_before(text, anchor, addition):
    if addition.strip() not in text:
        assert anchor in text, anchor
        text = text.replace(anchor, addition + "\n\n" + anchor, 1)
    return text


def replace_paragraph(text, prefix, replacement):
    parts = text.split("\n\n")
    matches = [i for i, p in enumerate(parts) if p.startswith(prefix)]
    assert len(matches) == 1, (prefix, matches)
    parts[matches[0]] = replacement
    return "\n\n".join(parts)


def save_ledger(root, name, additions, updates):
    path = root / "sources" / (name + "_sources.json")
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    rows = data if isinstance(data, list) else data["sources"]
    key = "id" if name == "heterogeneity" else "key"
    for row in rows:
        if row[key] in updates:
            row.update(updates[row[key]])
    existing = {r[key] for r in rows}
    rows.extend(r for r in additions if r[key] not in existing)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.with_suffix(".csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, type=Path)
    root = p.parse_args().root
    report = root / "reports" / "01_heterogeneity_review.md"
    text = report.read_text(encoding="utf-8-sig").replace("`SOURCE_LEDGER.md`", "`../sources/heterogeneity_source_ledger.md`")
    qian = r"""**A particularly close sparse competitor is [Qian, Chakraborty, Maiti and Cheung (2021)](https://doi.org/10.5705/ss.202018.0451).** Their semiparametric model has flexible prognosis and a linear treatment contrast, with an initial global null that all interaction slopes vanish. They screen marginal interactions and proceed through forward tests, using an adaptive $m$-out-of-$n$ bootstrap to address nonregularity. The paper develops a randomized-trial procedure and a doubly robust observational extension under stated nuisance conditions. Its sequential terminology concerns adding covariates, not optional-stopping-valid monitoring. This is direct prior art for sparse interaction testing. Our joint quadratic statistic and proposed subset-mixture evidence differ in how directions are combined; arbitrary same-sample forest fitting is not automatically covered by their theorem. [Accepted manuscript](https://www3.stat.sinica.edu.tw/preprint/SS-2018-0451_Preprint.pdf)."""
    text = insert_before(text, "Our recommendation is to retain", qian)
    selective = r"""**Selection changes the inferential problem.** [Zhao, Small and Ertefaie (2022)](https://doi.org/10.1111/rssb.12483) combine Robinson residualization, machine-learned nuisance estimates and lasso selection of an interpretable effect-modification model. They derive conditional selective inference under nuisance-rate assumptions. This addresses inference after choosing modifiers, whereas a prespecified global test asks whether an entire fixed collection is inactive. Neither reporting selected ordinary intervals nor inspecting a few posterior coefficients reproduces selection-adjusted inference. Their simulations also caution that observational nuisance estimation can challenge selective error control. This is a useful second-stage comparator if the project progresses from detecting heterogeneity to naming its drivers. [Accepted author manuscript](https://arxiv.org/abs/1705.08020); [indexed full-text discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC9491375/)."""
    text = insert_before(text, "## 8. Variance of CATE", selective)
    text = text.replace("develop a kernel-machine score test for markers predictive of treatment differences.", "develop a kernel-machine score test for markers predictive of treatment differences in randomized-trial settings.")
    text = text.replace("Kernel choice and effect scale must be matched before interpreting comparisons.", "Kernel choice, trial design and effect scale must be matched before interpreting comparisons. Applying the construction to observational data with estimated, varying propensity scores requires a separately justified adaptation.")
    text = text.replace("develop conditional qualitative-effect testing. These belong", "test whether additional covariates improve optimal treatment decisions conditional on an existing covariate set. That incremental decision target is more specific than merely asking whether a CATE takes both signs. These belong")
    text = replace_paragraph(text, "The [WATCH repository]", "The [WATCH workflow, Sechidis et al. (2025)](https://doi.org/10.1002/pst.2463), and its [official repository](https://github.com/Novartis/WATCH) organize planning, data preparation, exploration and multidisciplinary assessment. The final journal article was first online on 26 December 2024 and appears in the 2025 journal record. WATCH is explicitly exploratory and is broader than one global test; a local DR/conditional-inference reconstruction is not the whole workflow. Its contribution here is procedural: formulate scientific modifiers and effect scales before exploration, then separate evidence of variation from interpretation of discovered groups.")
    row = "| Qian et al. sequential interactions | Global and additional linear effect modifiers | Adaptive m-out-of-n bootstrap | Sparse semiparametric interaction testing | Forward selection and nuisance assumptions | Essential |"
    if row not in text:
        text = text.replace("| Shen–Cai kernel score", row + "\n| Shen–Cai kernel score", 1)
    text = text.replace("| Kernel and scale alignment | Essential |", "| Kernel, design and scale alignment | Essential |")
    text = insert_before(text, "Semi-BART is an additional", "Qian et al.'s testing construction further limits novelty claims about detecting sparse linear effect modifiers with a flexible prognostic component. The new contribution must be the particular evidence target, prior, computation or theorem, rather than the architecture alone.")
    report.write_text(text, encoding="utf-8")

    report = root / "reports" / "02_bayes_factor_literature_review.md"
    text = report.read_text(encoding="utf-8-sig")
    text = replace_paragraph(text, "**Datta,", "**Datta, Guha, Shudde and Johnson (2025)** combine statistic-based Bayes factors with nonlocal moment priors for z, t, chi-square and F statistics. The final four-author article supersedes the older three-author preprint and is direct prior art for a nonlocal transformation of a heterogeneity statistic. Its formulas concern specified sampling laws; an estimated forest score still needs justification. The [2025 article](https://doi.org/10.1016/j.spl.2024.110330) has a [corrigendum](https://doi.org/10.1016/j.spl.2025.110446) correcting the two-sided t formula in Table 1 and the supplement. Any implementation of that branch must use the correction. The Gaussian-score formula used in our pilot is derived separately.")
    bff = "**Bayes factor functions make effect-scale sensitivity explicit.** [Johnson, Pramanik and Shudde (2023)](https://doi.org/10.1073/pnas.2217331120) present evidence as a function of a scientifically interpretable effect-size index. [Datta et al., *On Bayes Factor Functions* (2025)](https://doi.org/10.1214/25-BA1562), extend this direction with inverse-moment alternatives for common statistics; the [author manuscript](https://arxiv.org/abs/2506.16674) explains the distinction from local priors. For this project, plotting evidence across plausible heterogeneity scales is more informative than reporting one unexplained default. A curve is a sensitivity display: choosing its largest value after observing the data does not turn it into the BF of a prespecified mixture. The latter integrates over a declared prior on scale. These are established precedents for any proposed evidence-versus-effect-size visualization."
    text = insert_before(text, "## 8. Closest model-building", bff)
    goeman = r"""**The variance-component connection is also established.** [Goeman, van de Geer, de Kort and van Houwelingen (2004)](https://doi.org/10.1093/bioinformatics/btg382) developed a global association test for groups of genes, not a causal heterogeneity test. The mathematical bridge is our synthesis: integrating a Gaussian modifier prior adds $AVA^\top$ to the outcome covariance. A variance-component score and a Gaussian integrated BF can therefore probe closely related covariance departures, while using different evidence scales. This is another reason that a Gaussian integrated alternative alone is insufficient novelty."""
    text = insert_before(text, "## 9. Theory:", goeman)
    report.write_text(text, encoding="utf-8")

    report = root / "reports" / "03_standalone_research_proposal.md"
    text = report.read_text(encoding="utf-8-sig")
    text = insert_before(text, "This proposal recommends proceeding", "The frequentist positioning must also include [Qian et al. (2021)](https://doi.org/10.5705/ss.202018.0451), whose flexible-prognosis/linear-treatment model tests global and forward-selected effect modifiers. A full-model Bayesian evidence calculation is a different inferential target, but the sparse semiparametric interaction-testing problem is established. The companion field review also covers selection-adjusted inference for discovered modifiers.")
    coding = r"""**Treatment coding is part of the model specification.** This proposal fixes $W=Z-e$ and a shared prior on the corresponding prognostic function. Rewriting $\mu(X)+Z\tau(X)$ as $m(X)+(Z-e)\tau(X)$ shifts prognosis by $e(X)\tau(X)$. Keeping an independent forest prior unchanged after that transformation generally changes the joint prior and hence the Bayes factor. The existing API also contains uncentered or adaptive coding and propensity updates. Those options are not automatically equivalent to the present experiment; any extension must define its induced priors and evidence target explicitly."""
    text = insert_before(text, "The null is conditional on", coding)
    energy = r"""**Separate sparsity from total prior effect variation.** With standardized approximately orthogonal modifier columns and independent slab variance $s^2$, a dense Gaussian prior has $E\|\delta\|^2=qs^2$. Under independent inclusion probability $\rho$ conditioned on a nonempty subset, it is $s^2q\rho/[1-(1-\rho)^q]$. At $q=6$, $\rho=1/6$, $s=0.5$, these are 1.5 and approximately 0.376. Equal slab scales therefore compare complete prior specifications, not sparsity alone. Add a prespecified energy-matched sensitivity comparison before attributing performance differences to sparse support. For correlated columns, match $E[\delta^\top\Sigma_h\delta]$ on a declared target covariate distribution."""
    text = insert_before(text, "Shared priors require special care", energy)
    text = insert_before(text, "## 4. Exact conditional integration", "The Gaussian covariance alternative below is related to the variance-component/global-test lineage, including [Goeman et al. (2004)](https://doi.org/10.1093/bioinformatics/btg382). The contribution sought is reliable nuisance integration and sparse evidence with matched theory, rather than a new Gaussian integration identity.")
    report.write_text(text, encoding="utf-8")

    def hrow(id, title, authors, year, doi, url, access, role):
        return dict(id=id, title=title, authors=authors, year=year, doi=doi, url=url, access_level=access, role=role, notes="Added after independent cross-review", verified_on="2026-10-16")
    def brow(key, title, authors, year, doi, url, access, role):
        return dict(key=key, title=title, authors=authors, year=year, doi=doi, url=url, access_level=access, role=role, accessed="2026-10-16")
    save_ledger(root, "heterogeneity", [
        hrow("qian2021", "A Sequential Significance Test for Treatment by Covariate Interactions", "Min Qian; Bibhas Chakraborty; Raju Maiti; Ying Kuen Cheung", 2021, "10.5705/ss.202018.0451", "https://www3.stat.sinica.edu.tw/preprint/SS-2018-0451_Preprint.pdf", "Primary accepted PDF methods/theorem sections inspected; final citation checked", "Direct sparse semiparametric interaction-testing competitor"),
        hrow("zhao2022selective", "Selective inference for effect modification via the lasso", "Qingyuan Zhao; Dylan S. Small; Ashkan Ertefaie", 2022, "10.1111/rssb.12483", "https://arxiv.org/abs/1705.08020", "Accepted author abstract and indexed PMC discussion checked; direct PMC blocked", "Conditional inference after modifier selection")
    ], {"watch2024": dict(year=2025, doi="10.1002/pst.2463", url="https://onlinelibrary.wiley.com/doi/abs/10.1002/pst.2463", notes="Final 2025 article; first online 26 December 2024; exploratory workflow", access_level="Final publisher abstract/metadata and official repository checked")})
    # Resolve the older Datta key without assuming its spelling.
    bdata = json.loads((root / "sources/bayes_sources.json").read_text(encoding="utf-8-sig"))
    brows = bdata if isinstance(bdata, list) else bdata["sources"]
    datta_key = next(r["key"] for r in brows if "Datta" in r.get("authors", "") and ("2310.16213" in r.get("url", "") or r.get("doi") == "10.1016/j.spl.2024.110330"))
    datta_authors = "Saptati Datta; Riana Guha; Rachael Shudde; Valen E. Johnson"
    save_ledger(root, "bayes", [
        brow("datta_corrigendum_2025", "Corrigendum to Bayes factors functions based on test statistics and non-local moment prior densities", datta_authors, 2025, "10.1016/j.spl.2025.110446", "https://www.sciencedirect.com/science/article/pii/S0167715225000914", "Primary corrigendum indexed text inspected by independent reviewer; root endpoint 403", "Correction to two-sided t formula; not used by Gaussian pilot"),
        brow("johnson_pramanik_shudde_2023", "Bayes factor functions for reporting outcomes of hypothesis tests", "Valen E. Johnson; Sandipan Pramanik; Rachael Shudde", 2023, "10.1073/pnas.2217331120", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9974512/", "Primary full text inspected by independent reviewer; root direct PMC challenge", "Evidence as function of effect-size index"),
        brow("datta_bff_2025", "On Bayes Factor Functions", datta_authors, 2025, "10.1214/25-BA1562", "https://arxiv.org/abs/2506.16674", "Primary author manuscript abstract and author publication list; DOI metadata fetched separately", "Inverse-moment alternatives and effect-scale evidence curves"),
        brow("goeman_2004", "A global test for groups of genes: testing association with a clinical outcome", "Jelle J. Goeman; Sara A. van de Geer; Floor de Kort; Hans C. van Houwelingen", 2004, "10.1093/bioinformatics/btg382", "https://academic.oup.com/bioinformatics/article-abstract/20/1/93/229017", "Primary publisher abstract/metadata checked by independent reviewer", "Adjacent variance-component association-test lineage")
    ], {datta_key: dict(title="Bayes factors functions based on test statistics and non-local moment prior densities", authors=datta_authors, year=2025, doi="10.1016/j.spl.2024.110330", url="https://doi.org/10.1016/j.spl.2024.110330", access_level="Final publisher metadata and preprint checked; separate corrigendum recorded")})
    for name in ["heterogeneity", "bayes"]:
        ledger = root / "sources" / (name + "_source_ledger.md")
        text = ledger.read_text(encoding="utf-8-sig")
        marker = "## Final independent-audit additions"
        if marker not in text:
            text += "\n\n" + marker + "\n\nThe machine-readable `" + name + "_sources.csv` and `.json` contain the final records and access limitations. Cross-review added Qian's sparse semiparametric test and Zhao's selective inference to the field review; updated WATCH to its final 2025 record; updated Datta's four-author paper and corrigendum; and added Bayes factor functions and the Goeman variance-component lineage to the Bayesian review. Earlier preprint labels above are version provenance, superseded by the final records where indicated. See the two preserved review memos in this folder.\n"
        ledger.write_text(text, encoding="utf-8")
    raw = json.loads((root / "grok/writing_raw.json").read_text(encoding="utf-8-sig"))
    output = raw.get("response", raw.get("text", raw.get("result", "")))
    if isinstance(output, str) and output:
        (root / "grok/WRITING_AIDE_UNVERIFIED.md").write_text("# Grok writing aide — unverified external draft\n\nRetained for provenance; final reports use independently checked wording and actual pilot specifications.\n\n" + output, encoding="utf-8")
    audit = root / "grok/GROK_USE_AND_CORRECTIONS.md"
    text = audit.read_text(encoding="utf-8-sig")
    if "## Writing-aide pass" not in text:
        text += "\n\n## Writing-aide pass\n\nA second read-only Grok CLI run used a supplied, corrected technical brief to suggest accessible explanations, critique ambiguous wording and propose organization. Its prognosis-versus-treatment example and prior-odds framing informed the introductory guide. Its warnings about finite Monte Carlo averages and finite modifier bases were retained. The brief did not include every prior parameter; criticisms of unspecified leaf variance, inclusion probability and covariance were checked against the actual code, which specifies them. Its assertion that the listed scales do not apply to the dense prior was not adopted: the pilot uses the same stated scale grid for both. The advice to leave results empty applied before simulation completion and does not govern the final empirical report. Raw prompts and responses are preserved; Grok is an aide, not a cited scientific authority.\n"
    audit.write_text(text, encoding="utf-8")
    print("Applied final audit corrections and source-ledger updates.")


if __name__ == "__main__":
    main()
