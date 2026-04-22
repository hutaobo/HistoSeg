$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
git -C $repoRoot config core.hooksPath .githooks
$hooksPath = git -C $repoRoot config --get core.hooksPath
Write-Output "Configured core.hooksPath=$hooksPath"
