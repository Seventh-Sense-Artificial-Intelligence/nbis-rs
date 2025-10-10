# PowerShell script to patch maturin wheel
# Alternative to patch_maturin_wheel.bat with better error handling

$ErrorActionPreference = "Stop"

Write-Host "Patching maturin wheel..." -ForegroundColor Green

$zlib1Path = "C:\msys64\mingw64\bin\zlib1.dll"
if (-not (Test-Path $zlib1Path)) {
    Write-Error "zlib1.dll not found at: $zlib1Path"
    exit 1
}

# Get the script's directory and look for dist directory there
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$distDir = Join-Path $scriptDir "dist"

# Change to dist directory
if (-not (Test-Path $distDir)) {
    Write-Error "dist directory not found at: $distDir"
    exit 1
}

Set-Location $distDir

# Find the wheel file
$wheelFiles = Get-ChildItem -Filter "*.whl"
if ($wheelFiles.Count -eq 0) {
    Write-Error "No wheel file found in dist directory"
    exit 1
}

$wheelFile = $wheelFiles[0]
Write-Host "Found wheel file: $($wheelFile.Name)" -ForegroundColor Cyan

# Get absolute path of wheel file
$wheelPath = $wheelFile.FullName

# Create temporary directory
$wheelUnzipDir = Join-Path $env:TEMP "wheel_patch_$(Get-Random)"
New-Item -ItemType Directory -Path $wheelUnzipDir | Out-Null

Write-Host "Created temporary directory: $wheelUnzipDir" -ForegroundColor Green

try {
    # Load .NET compression assembly
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    
    # Extract wheel file (wheel files are ZIP files)
    Write-Host "Extracting wheel file..." -ForegroundColor Cyan
    [System.IO.Compression.ZipFile]::ExtractToDirectory($wheelPath, $wheelUnzipDir)

    # Find and patch METADATA file
    Write-Host "🚫 Removing License-File: from METADATA..." -ForegroundColor Yellow
    $metadataFile = Get-ChildItem -Path $wheelUnzipDir -Name "METADATA" -Recurse | Select-Object -First 1
    
    if (-not $metadataFile) {
        Write-Error "METADATA file not found in wheel"
        exit 1
    }
    
    $metadataPath = Join-Path $wheelUnzipDir $metadataFile
    
    # Filter out License-File lines
    $content = Get-Content $metadataPath | Where-Object { $_ -notmatch '^License-File:' }
    $utf8NoBom = New-Object System.Text.UTF8Encoding $false
    [System.IO.File]::WriteAllLines($metadataPath, $content, $utf8NoBom)

    # Rename libuniffi_nbis.* → libnbis.*
    Write-Host "Renaming library files..." -ForegroundColor Cyan
    $nbisDir = Join-Path $wheelUnzipDir "nbis"

    Write-Host "Checking for nbis directory: $nbisDir" -ForegroundColor Cyan
    
    if (Test-Path $nbisDir) {
        Get-ChildItem -Path $nbisDir -Name "uniffi_nbis.*" | ForEach-Object {
            $oldPath = Join-Path $nbisDir $_
            # $newName = $_ -replace '^libuniffi_nbis', 'libnbis'
            $newName = $_ -replace '^uniffi_nbis', 'nbis'
            $newPath = Join-Path $nbisDir $newName

            Write-Host "Found library file: $_" -ForegroundColor Cyan
            Write-Host "Renaming: $oldPath → $newPath" -ForegroundColor White

            Write-Host "Renaming: $_ → $newName" -ForegroundColor White
            Move-Item $oldPath $newPath
        }
    } else {
        Write-Host "Warning: nbis directory not found, skipping library file renaming" -ForegroundColor Yellow
    }

    #copy zlib1.dll to the wheel directory
    Write-Host "Copying zlib1.dll to wheel directory..." -ForegroundColor Cyan
    $zlibDestPath = Join-Path $wheelUnzipDir "nbis"
    if (-not (Test-Path $zlibDestPath)) {
        New-Item -ItemType Directory -Path $zlibDestPath | Out-Null
    }
    Copy-Item -Path $zlib1Path -Destination $zlibDestPath -Force    

    # Remove the original wheel file
    Write-Host "Removing original wheel file..." -ForegroundColor Cyan
    Remove-Item $wheelPath

    #### Repack the wheel
    ####Write-Host "Repacking wheel file..." -ForegroundColor Cyan
    
    
    #### Create the new wheel file using .NET compression
    ####[System.IO.Compression.ZipFile]::CreateFromDirectory($wheelUnzipDir, $wheelPath)

    # Repack the wheel using Python (ensures forward slashes inside ZIP)
    Write-Host "Repacking wheel using Python..." -ForegroundColor Cyan

    $tempPy = Join-Path $env:TEMP "repack_wheel.py"

    # Write Python script line by line (no PowerShell indentation preserved)
    @(
    'import os, zipfile, sys'
    'root = sys.argv[1]'
    'wheel_path = sys.argv[2]'
    'print(f"Repacking wheel from {root} -> {wheel_path}")'
    'with zipfile.ZipFile(wheel_path, "w", zipfile.ZIP_DEFLATED) as zf:'
    '    for dirpath, _, filenames in os.walk(root):'
    '        for filename in filenames:'
    '            full_path = os.path.join(dirpath, filename)'
    '            rel_path = os.path.relpath(full_path, root).replace("\\\\", "/")'
    '            zf.write(full_path, rel_path)'
    'print("Repacking complete with forward slashes.")'
    ) | Set-Content -Path $tempPy -Encoding UTF8

    # Run Python with safe arguments
    python "$tempPy" "$wheelUnzipDir" "$wheelPath"

    if ($LASTEXITCODE -ne 0) {
        Write-Error "Python repacking failed"
        exit 1
    }

    Remove-Item $tempPy -Force
    Write-Host "✅ Wheel repacked successfully." -ForegroundColor Green

    Write-Host "✅ Patched and rebuilt wheel: dist\$($wheelFile.Name)" -ForegroundColor Green

} catch {
    Write-Error "Error during wheel patching: $($_.Exception.Message)"
    exit 1
} finally {
    # Return to original directory and cleanup
    Set-Location $scriptDir
    
    # if (Test-Path $wheelUnzipDir) {
    #     Remove-Item -Recurse -Force $wheelUnzipDir
    # }
}
