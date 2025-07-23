#!/bin/bash
set -e

# Set path to jniLibs target dir
JNILIBS_DIR="bindings/android/app/src/main/jniLibs"

# Android targets you want to build for
TARGETS=(
  "aarch64-linux-android"    # arm64-v8a
  "armv7-linux-androideabi"  # armeabi-v7a
)

# Build each target
for TARGET in "${TARGETS[@]}"; do
  echo "📦 Building for $TARGET..."
  cargo ndk -t "$TARGET" -o ./target/android build --release
done

# Clear jniLibs output
echo "🧹 Cleaning existing jniLibs..."
rm -f "$JNILIBS_DIR/arm64-v8a/"*.so
rm -f "$JNILIBS_DIR/armeabi-v7a/"*.so

# Copy .so files into correct ABI folders
for TARGET in "${TARGETS[@]}"; do
  case "$TARGET" in
    "aarch64-linux-android")
      ABI="arm64-v8a"
      ;;
    "armv7-linux-androideabi")
      ABI="armeabi-v7a"
      ;;
    *)
      echo "❌ Unknown target: $TARGET"
      exit 1
      ;;
  esac

  echo "📁 Copying .so for $ABI..."
  mkdir -p "$JNILIBS_DIR/$ABI"

  # Adjust the library name as needed
  LIB_NAME="libnbis.so"
  cp "target/android/$ABI/$LIB_NAME" "$JNILIBS_DIR/$ABI/"
done

echo "✅ Done. .so files copied to $JNILIBS_DIR"

# Detect platform-specific dynamic library extension
case "$(uname)" in
  Darwin)   LIB_EXT="dylib" ;;
  Linux)    LIB_EXT="so" ;;
  *)        echo "❌ Unsupported platform: $(uname)" && exit 1 ;;
esac

cargo build --release
cargo run --bin uniffi-bindgen generate \
          --library "target/release/libnbis.${LIB_EXT}" \
          --language kotlin \
          --out-dir bindings/android/app/src/main/java/com/sensecrypt/sdk/nbis

# Move to sensible folder without uniffi opinionation
mv bindings/android/app/src/main/java/com/sensecrypt/sdk/nbis/uniffi/nbis/nbis.kt \
    bindings/android/app/src/main/java/com/sensecrypt/sdk/nbis.kt
rm -rf bindings/android/app/src/main/java/com/sensecrypt/sdk/nbis


# Root directory of generated Kotlin files
TARGET_DIR="bindings/android/app/src/main/java/com/sensecrypt/sdk"

# Old and new package names
OLD_PACKAGE="package uniffi.nbis"
NEW_PACKAGE="package ai.seventhsense.sdk.nbis"

# Recursively find all .kt files and replace the package line
echo "🔍 Replacing package '$OLD_PACKAGE' with '$NEW_PACKAGE' in $TARGET_DIR"

# Use sed with platform detection
if [[ "$OSTYPE" == "darwin"* ]]; then
  # macOS sed requires an empty string after -i
  find "$TARGET_DIR" -name "*.kt" -exec sed -i '' "s|^$OLD_PACKAGE|$NEW_PACKAGE|" {} +
else
  # Linux sed
  find "$TARGET_DIR" -name "*.kt" -exec sed -i "s|^$OLD_PACKAGE|$NEW_PACKAGE|" {} +
fi

# Build the aar file
echo "📦 Building Android AAR..."
cd bindings/android
./gradlew :app:assembleRelease
if [ $? -ne 0 ]; then
  echo "❌ Failed to build Android AAR"
  exit 1
fi
# Rename the AAR file
AAR_FILE=$(find app/build/outputs/aar -name "*.aar" | head -n 1)
if [ -z "$AAR_FILE" ]; then
  echo "❌ No AAR file found"
  exit 1
fi
mv "$AAR_FILE" "../../dist/nbis-sdk.aar"
