$projectRoot = Split-Path -Parent $PSScriptRoot
& wsl.exe -d Ubuntu-24.04 --cd $projectRoot -- bash --rcfile scripts/moose_bashrc -i
exit $LASTEXITCODE
