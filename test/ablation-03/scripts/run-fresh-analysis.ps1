param(
  [string]$CacheRoot = 'D:\cache\ccs\_ablation-03-fresh',
  [int]$TotalCores = 16,
  [int]$Workers = 2,
  [int]$MemoryGB = 48,
  [switch]$SkipBenchmark
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
$rscript = 'C:\R\R-4.3.1\bin\Rscript.exe'
if (-not (Test-Path -LiteralPath $rscript)) {
  throw "Windows R 4.3.1 not found: $rscript"
}
if ($TotalCores -lt 1 -or $Workers -lt 1 -or $Workers -gt $TotalCores) {
  throw 'TotalCores must be positive and Workers must be between 1 and TotalCores.'
}
if ($MemoryGB -lt 1) { throw 'MemoryGB must be positive.' }

$CacheRoot = [IO.Path]::GetFullPath($CacheRoot)
New-Item -ItemType Directory -Path $CacheRoot -Force | Out-Null
$lock = Join-Path $CacheRoot '.workflow-lock'
if (Test-Path -LiteralPath $lock) {
  throw "Cache lock exists: $lock. Confirm no R process is active, then remove the lock manually before retrying."
}

$env:LANG = 'Chinese_China.utf8'
$env:LC_ALL = 'Chinese_China.utf8'
$env:LC_CTYPE = 'Chinese_China.utf8'
$env:CCS_ABLATION_CACHE_ROOT = $CacheRoot
$env:CCS_ABLATION_CORES = [string]$TotalCores
$env:CCS_ABLATION_WORKERS = [string]$Workers
$env:CCS_ABLATION_MEMORY_GB = [string]$MemoryGB

Set-Location $repoRoot
$stages = @(
  'test/ablation-03/01.00.00. 数据准备.R',
  'test/ablation-03/02.00.00. 表示输入准备.R',
  'test/ablation-03/03.00.00. 生物输入准备.R'
)
if (-not $SkipBenchmark) {
  $stages += 'test/ablation-03/scripts/benchmark-learning-curve.R'
}
$stages += @(
  'test/ablation-03/05.00.00. 表示分析.R',
  'test/ablation-03/06.00.00. 生物锚点分析.R',
  'test/ablation-03/07.00.00. 结构复现分析.R'
)

Write-Host "Fresh ablation-03 cache: $CacheRoot"
Write-Host "CPU budget: $Workers worker(s) x up to $([math]::Floor($TotalCores / $Workers)) XGBoost thread(s); memory budget: ${MemoryGB} GB"
foreach ($stage in $stages) {
  Write-Host "`n>>> $stage"
  & $rscript --vanilla $stage
  if ($LASTEXITCODE -ne 0) { throw "Stage failed with exit code $LASTEXITCODE`: $stage" }
}

Write-Host "`nFresh ablation-03 analysis completed. Cache root: $CacheRoot"
