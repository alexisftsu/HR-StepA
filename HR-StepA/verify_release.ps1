param(
  [string]$SumsPath = ".\out\SHA256SUMS.txt",
  [string[]]$Assets = @(
    ".\out\HR_Final_pkg_FINAL_BEST_beta0.08_D1100.zip",
    ".\out\HR_Final_inputs_beta0.08_D1100.zip",
    ".\lemma_bs_insertion_filled.pdf",
    ".\out\BEST_summary.txt",
    ".\out\summary_sweep_beta0.08.csv",
    ".\out\summary_sweep_beta0.09.csv",
    ".\out\summary_sweep_beta0.1.csv",
    ".\out\summary_sweep_beta0.11.csv",
    ".\out\summary_sweep_beta0.12.csv"
  )
)
$hashMap = @{}
Get-Content $SumsPath | ForEach-Object {
  if ($_ -match '^\s*([0-9A-Fa-f]{64})\s+(.+)$') {
    $hashMap[$Matches[2]] = $Matches[1].ToUpper()
  }
}
$ok=$true
foreach ($p in $Assets) {
  $name = Split-Path $p -Leaf
  if (-not (Test-Path $p)) { Write-Host "MISSING: $name" -ForegroundColor Red; $ok=$false; continue }
  $h = (Get-FileHash $p -Algorithm SHA256).Hash.ToUpper()
  $ref = $hashMap[$name]
  if ($ref -and $ref -eq $h) { Write-Host "OK  $name" -ForegroundColor Green }
  else { Write-Host "FAIL $name`n  got: $h`n  exp: $ref" -ForegroundColor Red; $ok=$false }
}
if (-not $ok) { exit 1 }
