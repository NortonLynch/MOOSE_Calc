$projectRoot = Split-Path -Parent $PSScriptRoot
& wsl.exe -d Ubuntu-24.04 --cd $projectRoot -- bash scripts/verify_environment.sh
exit $LASTEXITCODE
