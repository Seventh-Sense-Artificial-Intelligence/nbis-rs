#!/usr/bin/env bash
set -euo pipefail

# On Debian/Ubuntu, ensure python3-venv is installed so `python3 -m venv` will work
if [ -f /etc/os-release ] && grep -qiE 'ubuntu|debian' /etc/os-release; then
  echo "Detected Debian/Ubuntu; installing python3-venv..."
  sudo apt-get update
  sudo apt-get install -y python3-venv
fi

python3 -m venv .venv
source .venv/bin/activate
pip install maturin
maturin build --release --manifest-path Cargo.toml --out dist
./patch_maturin_wheel.sh