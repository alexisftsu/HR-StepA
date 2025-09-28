param(
  [Parameter(Mandatory=$true)][string]$Tag,
  [string]$Repo = "alexisftsu/HR-StepA",
  [switch]$RunVerify  # dispara el workflow de checksums
)

$ErrorActionPreference = 'Stop'

# Ir a HR-StepA dentro del repo
$root = git rev-parse --show-toplevel | % Trim
Set-Location "$root/HR-StepA"

# 1) Validar manifest
& "$root/tools/Test-BestManifest.ps1" -Path "$PWD/BEST_manifest.json"

# 2) Cargar lista de assets desde el manifest (quitar 'out/' si viene)
$m = Get-Content -Raw BEST_manifest.json | ConvertFrom-Json
$assets = @(
  $m.assets.best_zip,
  $m.assets.inputs_zip,
  $m.assets.lemma_pdf,
  'BEST_summary.txt'
) | % { ($_ -replace '^out/','') } | Where-Object { $_ -and (Test-Path $_) }

if (-not $assets -or $assets.Count -lt 3) {
  throw "Faltan assets requeridos en el directorio actual: $PWD"
}

# 3) Crear release si no existe
try { gh release view $Tag --repo $Repo | Out-Null }
catch { gh release create $Tag --repo $Repo -t $Tag -n "HR release $Tag" }

# 4) Subir assets (idempotente)
foreach ($a in $assets) {
  gh release upload $Tag $a --repo $Repo --clobber
}

# 5) Subir SHA256SUMS.txt si existe (si no, lo generará el workflow)
if (Test-Path "$PWD/SHA256SUMS.txt") {
  gh release upload $Tag "$PWD/SHA256SUMS.txt" --repo $Repo --clobber
}

# 6) Disparar verificación (Linux genera sums si faltan; Windows verifica)
if ($RunVerify) {
  gh workflow run ".github/workflows/release-checksums.yml" -f tag=$Tag --ref main --repo $Repo | Out-Null
  "Workflow lanzado para $Tag."
}

"Release $Tag actualizado en $Repo."
