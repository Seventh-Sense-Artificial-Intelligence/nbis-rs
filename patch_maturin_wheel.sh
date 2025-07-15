#!/bin/bash
set -euo pipefail

pushd dist > /dev/null

# Find the wheel file and get its absolute path
wheel_file=$(ls *.whl | head -n 1)
wheel_path=$(realpath "$wheel_file")
wheel_unzip_dir=$(mktemp -d)

# Unzip wheel into temp dir
unzip "$wheel_file" -d "$wheel_unzip_dir"

echo "🚫 Removing License-File: from METADATA..."
metadata_file=$(find "$wheel_unzip_dir" -name METADATA)
tmpfile=$(mktemp)
grep -v '^License-File:' "$metadata_file" > "$tmpfile"
mv "$tmpfile" "$metadata_file"

# Rename libuniffi_nbis.* → libnbis.*
for f in "$wheel_unzip_dir"/nbis/libuniffi_nbis.*; do
    mv "$f" "${f/libuniffi_nbis/libnbis}"
done

# Remove the original wheel
rm "$wheel_file"

# Repack the wheel correctly
cd "$wheel_unzip_dir"
shopt -s dotglob
zip -r "$wheel_path" * > /dev/null

popd > /dev/null
rm -rf "$wheel_unzip_dir"

echo "✅ Patched and rebuilt wheel: dist/$(basename "$wheel_path")"