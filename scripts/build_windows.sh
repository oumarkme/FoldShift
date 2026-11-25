#!/usr/bin/env bash
# Cross-compiles the Windows exes (x64 and ARM64) from Linux/WSL with cargo-xwin. Static C runtime: nothing to install.
# Usage: CARGO_TARGET_DIR=<local dir> scripts/build_windows.sh   (writes dist/FoldShift_v<version>_windows_{x64,arm64}.exe)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
mkdir -p dist
for pair in x86_64-pc-windows-msvc:x64 aarch64-pc-windows-msvc:arm64; do
  target=${pair%%:*}; tag=${pair##*:}
  rustup target add "$target" >/dev/null 2>&1 || true
  RUSTFLAGS="-C target-feature=+crt-static --remap-path-prefix=$root=/src --remap-path-prefix=$HOME=/home" \
    cargo xwin build --release --target "$target" -p psplot
  out="dist/FoldShift_v${version}_windows_${tag}.exe"
  cp "${CARGO_TARGET_DIR:-target}/$target/release/FoldShift.exe" "$out"
  (cd dist && sha256sum "$(basename "$out")" > "$(basename "$out").sha256")
  ls -la "$out"
done
