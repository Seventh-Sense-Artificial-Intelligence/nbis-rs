#!/usr/bin/env sh
set -e

##############################################################################
# 1. Build the Rust cdylib (same for all targets)
##############################################################################
cargo build --release

##############################################################################
# 2. Pick the right extension for the shared library -------------------------
#    Darwin → libxxx.dylib   |   Linux → libxxx.so
##############################################################################
case "$(uname -s)" in
  Darwin*) LIB_EXT="dylib" ;;
  Linux*)  LIB_EXT="so"    ;;
  *) echo "Unsupported OS: $(uname -s)" >&2; exit 1 ;;
esac

LIB_NAME="libnbis.${LIB_EXT}"
LIB_PATH="target/release/${LIB_NAME}"

##############################################################################
# 3. Run uniffi-bindgen to generate the Python stub
##############################################################################
cargo run --bin uniffi-bindgen generate \
  --library "${LIB_PATH}" \
  --language python \
  --out-dir bindings/python/nbispy

##############################################################################
# 4. Copy the compiled library next to the stub so Python can load it
##############################################################################
install -Dm755 "${LIB_PATH}" "bindings/python/nbispy/${LIB_NAME}"

echo "✔ Generated bindings/python/nbispy/{nbis.py,${LIB_NAME}}"