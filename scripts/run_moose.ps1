param(
    [ValidateRange(1, 20)][int]$Ranks = 1,
    [string]$InputFile,
    [string[]]$ExtraArgs = @()
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$arguments = @()
if ($InputFile) {
    $candidate = if ([System.IO.Path]::IsPathRooted($InputFile)) { $InputFile } else { Join-Path $projectRoot $InputFile }
    $resolved = (Resolve-Path -LiteralPath $candidate).Path
    if ($resolved -notmatch '^([A-Za-z]):[\\/](.*)$') { throw 'Input must be on a local Windows drive.' }
    $linuxInput = '/mnt/' + $Matches[1].ToLowerInvariant() + '/' + $Matches[2].Replace('\', '/')
    $arguments += @('-i', $linuxInput)
}
$arguments += $ExtraArgs
& wsl.exe -d Ubuntu-24.04 --cd $projectRoot -- bash scripts/run_moose.sh $Ranks @arguments
exit $LASTEXITCODE
