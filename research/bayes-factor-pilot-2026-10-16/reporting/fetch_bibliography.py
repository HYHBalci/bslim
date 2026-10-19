"""Fetch BibTeX from Crossref; never invent missing bibliographic fields.

Scans the delivered reports/source ledgers for DOI links or a doi CSV column.
Archives every raw response and records failures for manual follow-up.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import time
from urllib.parse import quote, unquote
from urllib.request import Request, urlopen


def collect(root):
    found = set()
    for folder in [root / "reports", root / "sources"]:
        for file in folder.glob("*"):
            if file.suffix.lower() not in {".md", ".csv", ".json"}:
                continue
            if file.name in {"crossref_responses.json", "bibliography_fetch_status.csv"}:
                continue
            text = file.read_text(encoding="utf-8-sig")
            if file.suffix.lower() == ".csv":
                for row in csv.DictReader(text.splitlines()):
                    for key, value in row.items():
                        if key and key.lower() == "doi" and value:
                            value = re.sub(r"^https?://(?:dx\.)?doi.org/", "", value.strip())
                            if value.startswith("10."):
                                found.add(value.rstrip("."))
                continue
            for doi in re.findall(r"https?://(?:dx\.)?doi\.org/(10\.[^\s<>\]\"},*]+)", text):
                doi = unquote(doi).rstrip(".,;`")
                if doi.endswith(")") and doi.count(")") > doi.count("("):
                    doi = doi[:-1]
                found.add(doi)
    return sorted(found, key=str.lower)


def fetch(doi):
    endpoints = ["https://api.crossref.org/works/" + quote(doi, safe="") + "/transform/application/x-bibtex",
                 "https://doi.org/" + quote(doi, safe="/")]
    if doi.lower().startswith("10.48550/"):
        endpoints = endpoints[1:]
    errors = []
    for url in endpoints:
        try:
            time.sleep(.7)  # Respect metadata-service rate limits.
            request = Request(url, headers={"Accept": "application/x-bibtex", "User-Agent": "ResearchBibliography/1.0 (literature review; metadata only)"})
            with urlopen(request, timeout=30) as response:
                text = response.read().decode("utf-8")
            if not text.lstrip().startswith("@"):
                raise ValueError("Response was not BibTeX")
            return {"doi": doi, "status": "fetched", "endpoint": url, "bibtex": text}
        except Exception as exc:
            errors.append(str(exc))
    return {"doi": doi, "status": "failed", "endpoint": endpoints[-1], "error": "; ".join(errors)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    a = parser.parse_args()
    dois = collect(a.root)
    sources = a.root / "sources"
    sources.mkdir(parents=True, exist_ok=True)
    cache_file = sources / "crossref_responses.json"
    cache = {r["doi"]: r for r in json.loads(cache_file.read_text(encoding="utf-8"))} if cache_file.exists() else {}
    records = [cache[d] if cache.get(d, {}).get("status") == "fetched" else fetch(d) for d in dois]
    (sources / "crossref_responses.json").write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    bibs = []
    for i, record in enumerate(records, 1):
        if record["status"] == "fetched":
            # Ensure a unique deterministic citation key without changing source metadata.
            key = f"verified_{i:03d}"
            bib = re.sub(r"(^\s*@\w+\s*\{)[^,]+,", lambda m: m.group(1) + key + ",", record["bibtex"], count=1)
            bibs.append("% Retrieved from " + record["endpoint"] + "\n" + bib.strip())
            record["citation_key"] = key
    (sources / "verified_references.bib").write_text("\n\n".join(bibs) + "\n", encoding="utf-8")
    with (sources / "bibliography_fetch_status.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["doi", "status", "citation_key", "endpoint", "error"], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps({"requested": len(dois), "fetched": len(bibs), "failed": [r["doi"] for r in records if r["status"] != "fetched"]}, indent=2))


if __name__ == "__main__":
    main()
