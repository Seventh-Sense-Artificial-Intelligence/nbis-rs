# PowerShell script to build Python wheel using maturin

param(
    [switch]$SkipVenvRecreate
)

$ErrorActionPreference = "Stop"

Write-Host "Building Python wheel on Windows (PowerShell)..." -ForegroundColor Green

# Get the absolute path dynamically
$projectPath = (Get-Location).Path
$env:RUSTFLAGS="-L $projectPath\ext\nfiq2_libs\x64\mingw\staticlib"

# Check if Python is available
try {
    $pythonVersion = python --version
    Write-Host "Found Python: $pythonVersion" -ForegroundColor Cyan
} catch {
    Write-Error "Python is not installed or not in PATH. Please install Python from https://www.python.org/downloads/"
    exit 1
}

# Check if Rust/Cargo is available
try {
    $cargoVersion = cargo --version
    Write-Host "Found Cargo: $cargoVersion" -ForegroundColor Cyan
} catch {
    Write-Error "Rust/Cargo is not installed. Please install Rust from https://rustup.rs/"
    exit 1
}

# Create or recreate virtual environment
if (Test-Path ".venv") {
    if (-not $SkipVenvRecreate) {
        Write-Host "Removing existing virtual environment..." -ForegroundColor Yellow
        Remove-Item -Recurse -Force ".venv"
    } else {
        Write-Host "Using existing virtual environment..." -ForegroundColor Cyan
    }
}

if (-not (Test-Path ".venv")) {
    Write-Host "Creating virtual environment..." -ForegroundColor Cyan
    python -m venv .venv
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Failed to create virtual environment. Make sure you have the full Python installation."
        exit 1
    }
}

# Activate virtual environment
Write-Host "Activating virtual environment..." -ForegroundColor Cyan
& ".venv\Scripts\Activate.ps1"

# Upgrade pip
Write-Host "Upgrading pip..." -ForegroundColor Cyan
python -m pip install --upgrade pip

# Install maturin
Write-Host "Installing maturin..." -ForegroundColor Cyan
pip install maturin
if ($LASTEXITCODE -ne 0) {
    Write-Error "Failed to install maturin"
    exit 1
}

# Create dist directory
if (-not (Test-Path "dist")) {
    New-Item -ItemType Directory -Name "dist" | Out-Null
}

# Build the wheel
Write-Host "Building wheel with maturin..." -ForegroundColor Cyan
maturin build --release --manifest-path Cargo.toml --out dist
if ($LASTEXITCODE -ne 0) {
    Write-Error "Maturin build failed"
    exit 1
}

# Run the patch script
Write-Host "Running wheel patch script..." -ForegroundColor Cyan
& ".\patch_maturin_wheel.ps1"
if ($LASTEXITCODE -ne 0) {
    Write-Error "Wheel patching failed"
    exit 1
}

Write-Host "✅ Build completed successfully!" -ForegroundColor Green
Write-Host "Wheel file is available in the dist/ directory" -ForegroundColor Cyan

# List the created wheel files
Write-Host "`nGenerated wheel files:" -ForegroundColor Yellow
Get-ChildItem -Path "dist" -Filter "*.whl" | ForEach-Object {
    Write-Host "  * $($_.Name)" -ForegroundColor White
}
