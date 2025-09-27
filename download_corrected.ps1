# Script de descarga corregido basado en análisis MD5
$fileListJson = Get-Content "file_list.json" | ConvertFrom-Json

Write-Host "Iniciando descarga de $($fileListJson.Count) archivos..." -ForegroundColor Yellow

# Resto del código de descarga aquí...
