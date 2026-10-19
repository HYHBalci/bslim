param(
  [string]$OutputDirectory = 'D:/heterogeneity-bayes-research-2026-10-16/pilot-results'
)
$ErrorActionPreference = 'Stop'
$repoRoot = (Get-Location).Path
$pilotRoot = Join-Path $repoRoot 'research/bayes-factor-pilot-2026-10-16'
$rBinary = 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe'
$nativeBinary = Join-Path $repoRoot '.research-library/fresh-2026-10-01/stochtree/libs/x64/stochtree.dll'
$allPaths = @(& rg --files R src include) + @(
  Join-Path $pilotRoot 'bayes_factor_sidecar.R'
  Join-Path $pilotRoot 'test_sidecar.R'
  Join-Path $pilotRoot 'run_pilot.R'
  Join-Path $pilotRoot 'summarize_pilot.R'
  Join-Path $pilotRoot 'analyze_results.R'
  Join-Path $pilotRoot 'run_diagnostics.R'
  Join-Path $pilotRoot 'publish_pilot_report.R'
  Join-Path $pilotRoot 'README.md'
  Join-Path $pilotRoot 'finalize_manifest.ps1'
  $rBinary
  $nativeBinary
  Join-Path $repoRoot '.research-library/fresh-2026-10-01/stochtree/R/stochtree.rdb'
  Join-Path $repoRoot '.research-library/fresh-2026-10-01/stochtree/R/stochtree.rdx'
  Join-Path $repoRoot '.research-library/fresh-2026-10-01/stochtree/DESCRIPTION'
  Join-Path $repoRoot '.research-library/fresh-2026-10-01/stochtree/NAMESPACE'
)
$hashes = foreach ($item in $allPaths) {
  if (Test-Path -LiteralPath $item -PathType Leaf) {
    $f = Get-Item -LiteralPath $item
    [pscustomobject]@{ path=$f.FullName; bytes=$f.Length; sha256=(Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash }
  }
}
$revision = & git rev-parse HEAD
$metadata = [ordered]@{
  captured_at=(Get-Date).ToString('o')
  project_root=$repoRoot
  git_head=$revision
  r_version='4.4.3'
  package_version='0.1.1'
  code_scope='Opt-in research sidecar; current R score wrappers sourced in an isolated environment; installed native binary reused'
  native_source_correspondence='Not established: hashes of current sources and installed DLL are distinct provenance records, not proof of a matching build'
  result_directory=$OutputDirectory
  source_and_binary_hashes=@($hashes)
}
$metadata | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $OutputDirectory 'provenance.json') -Encoding utf8
$hashes | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'source_and_binary_sha256.csv') -NoTypeInformation -Encoding utf8
Copy-Item -LiteralPath (Join-Path $pilotRoot 'README.md') -Destination (Join-Path $OutputDirectory 'METHODS_AND_REPRODUCTION.md') -Force
$testOutput = & $rBinary --vanilla (Join-Path $pilotRoot 'test_sidecar.R') 2>&1
$testOutput | Set-Content -LiteralPath (Join-Path $OutputDirectory 'mathematical_test_log.txt') -Encoding utf8
if ($LASTEXITCODE -ne 0) { throw 'Sidecar mathematical validation failed.' }
$reportPath = Join-Path (Split-Path $OutputDirectory -Parent) 'reports/04_pilot_results.md'
if (Test-Path -LiteralPath $reportPath) {
  Get-FileHash -LiteralPath $reportPath -Algorithm SHA256 | Select-Object Path,Hash | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $OutputDirectory 'report_sha256.json') -Encoding utf8
}
$artifactHashes = Get-ChildItem -LiteralPath $OutputDirectory -File | Where-Object { $_.Name -ne 'artifact_sha256.csv' } | ForEach-Object {
  [pscustomobject]@{name=$_.Name;bytes=$_.Length;sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash}
}
$artifactHashes | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'artifact_sha256.csv') -NoTypeInformation -Encoding utf8
Write-Output "Provenance, mathematical-test log, method copy, and artifact hashes saved to $OutputDirectory"
