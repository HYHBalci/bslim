"""Capture user-requested Grok research/editing assistance without shell encoding loss."""
import argparse
from pathlib import Path
import subprocess

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
p.add_argument("--prompt", required=True)
p.add_argument("--output", required=True)
p.add_argument("--writing-only", action="store_true")
a = p.parse_args()
command = [str(Path.home() / ".grok/bin/grok.exe"), "--cwd", str(a.root),
           "--permission-mode", "plan", "--no-subagents", "--max-turns", "3" if a.writing_only else "14",
           "--prompt-file", str(a.root / "grok" / a.prompt), "--output-format", "json",
           "--rules", "Return only your memo. Do not access local files, use shell tools, or change files. " +
           ("Use no tools; this is an editorial response to the supplied text." if a.writing_only else "Use only web search and web fetch for research.")]
if a.writing_only:
    command += ["--disable-web-search"]
with (a.root / "grok" / a.output).open("wb") as out, (a.root / "grok" / (a.output + ".stderr.txt")).open("wb") as err:
    result = subprocess.run(command, stdout=out, stderr=err)
raise SystemExit(result.returncode)
