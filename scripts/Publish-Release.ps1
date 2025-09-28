param(
  [Parameter(Mandatory=$true)][string]$Tag,
  [string]$Repo = "alexisftsu/HR-StepA",
  [switch]$RunVerify
)
$ErrorActionPreference = "Stop"

# Guard anti-tag de ejemplo
if ($Tag -match "beta0\.XX|DeltaYYYY") {
  Write-Warning "Tag de ejemplo detectado: $Tag. No haré descargas ni uploads."
  if ($RunVerify) {
    gh workflow run ".github/workflows/release-checksums.yml" -f tag=$Tag --ref main --repo $Repo | Out-Null
    "Workflow lanzado para $Tag (sin assets; Windows se saltará)."
  }
  return
}

# Ir a HR-StepA dentro del repo
$root = (git rev-parse --show-toplevel).Trim()
Set-Location "$root/HR-StepA"

# 1) Validar manifest
& "$root/tools/Test-BestManifest.ps1" -Path "$PWD/BEST_manifest.json" | Out-Host

# 2) Construir lista deseada + auto-descarga de faltantes
$m = Get-Content -Raw BEST_manifest.json | ConvertFrom-Json
$allWanted = @(
  $m.assets.best_zip,
  $m.assets.inputs_zip,
  $m.assets.lemma_pdf,
  "BEST_summary.txt"
) | ForEach-Object { ($_ -replace "^out/","") } | Where-Object { $_ }

$missing = $allWanted | Where-Object { -not (Test-Path $_) }
if ($missing) {
  foreach ($n in $missing) {
    gh release download $Tag --repo $Repo --skip-existing -p (Split-Path $n -Leaf) -D .
  }
}

$assets = $allWanted | Where-Object { Test-Path $_ }
$csvs = Get-ChildItem -Name "summary_sweep_beta*.csv" -ErrorAction SilentlyContinue
if ($csvs) { $assets += $csvs }

if ($assets.Where({$_ -match "\.zip$"}).Count -lt 2 -or -not ($assets -contains "lemma_bs_insertion_filled.pdf")) {
  throw "Faltan assets requeridos en el directorio actual: $PWD"
}

# 3) Crear release si no existe
$releaseExists = $true
try { gh release view $Tag --repo $Repo | Out-Null } catch { $releaseExists = $false }
if (-not $releaseExists) {
  gh release create $Tag --repo $Repo -t $Tag -n "HR release $Tag"
}

# 4) Subir assets (idempotente)
foreach ($a in $assets) {
  gh release upload $Tag $a --repo $Repo --clobber
}

# 5) Subir SHA256SUMS si existe
if (Test-Path "$PWD/SHA256SUMS.txt") {
  gh release upload $Tag "$PWD/SHA256SUMS.txt" --repo $Repo --clobber
}

# 6) Disparar verificación
if ($RunVerify) {
  gh workflow run ".github/workflows/release-checksums.yml" -f tag=$Tag --ref main --repo $Repo | Out-Null
  "Workflow lanzado para $Tag."
}
"Release $Tag actualizado en $Repo."
