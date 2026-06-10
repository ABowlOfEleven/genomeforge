#!/usr/bin/env bash
# Set GenomeForge's version everywhere it is declared, from a single command.
#
#   scripts/bump-version.sh 0.1.1
#
# Updates: Cargo workspace version, the WiX installer, the macOS Info.plist, and the
# AppStream metainfo. Then regenerates Cargo.lock's workspace entries. After running,
# update CHANGELOG.md, commit, and tag v<version>.
set -euo pipefail

v="${1:-}"
if [[ ! "$v" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "usage: scripts/bump-version.sh <x.y.z>   (e.g. 0.1.1)" >&2
  exit 1
fi

root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

# Cargo workspace version (drives every crate via [workspace.package]).
sed -i -E "s/^version = \"[0-9]+\.[0-9]+\.[0-9]+\"/version = \"$v\"/" Cargo.toml

# WiX MSI ProductVersion (four-part: x.y.z.0).
sed -i -E "s/Version=\"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+\"/Version=\"$v.0\"/" \
  installer/genomeforge.wxs

# macOS bundle (CFBundleVersion + CFBundleShortVersionString).
sed -i -E "s#<string>[0-9]+\.[0-9]+\.[0-9]+</string>#<string>$v</string>#g" \
  packaging/macos/Info.plist

# AppStream metainfo: the most recent <release> entry.
sed -i -E "0,/<release version=\"[0-9]+\.[0-9]+\.[0-9]+\"/s//<release version=\"$v\"/" \
  packaging/linux/io.github.abowlofeleven.GenomeForge.metainfo.xml

# Keep Cargo.lock's workspace member versions in sync (no dependency changes).
cargo update --workspace >/dev/null 2>&1 || true

echo "Version set to $v in:"
echo "  Cargo.toml, installer/genomeforge.wxs, packaging/macos/Info.plist,"
echo "  packaging/linux/*.metainfo.xml, Cargo.lock"
echo
echo "Next: edit CHANGELOG.md, commit, then:  git tag v$v && git push origin v$v"
