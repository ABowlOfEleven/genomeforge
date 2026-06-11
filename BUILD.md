# Building & packaging GenomeForge

## Run from source
```sh
cargo run -p gx-app                                   # debug
cargo run -p gx-app -- assets/sample_genome_23andme.txt   # auto-open a file
```

## Portable release executable
```sh
cargo build --release
# -> target/release/GenomeForge.exe   (single, self-contained, icon + metadata embedded)
```
The icon and version metadata are embedded automatically by `crates/gx-app/build.rs`
(via `winresource` + the Windows SDK `rc.exe`). If no resource compiler is found the
build still succeeds, just without the embedded icon.

## Windows MSI installer
Requires the free **WiX Toolset v5** (.NET tool; needs the .NET SDK):
```sh
dotnet tool install --global wix --version 5.0.2      # one-time
```
> Note: WiX v6/v7 require accepting a paid "Open Source Maintenance Fee" EULA, so stay on v5.

Then, from the repo root:
```sh
cargo build --release
mkdir -p target/wix
wix build installer/genomeforge.wxs -arch x64 -o target/wix/GenomeForge-0.1.0-x64.msi
# -> target/wix/GenomeForge-0.1.0-x64.msi
```
The installer is **per-user** (installs to `%LocalAppData%\GenomeForge`, no admin needed),
adds a Start-menu shortcut with the app icon, and supports clean upgrade/uninstall.
Install: double-click the `.msi`, or `msiexec /i target\wix\GenomeForge-0.1.0-x64.msi`.

## Regenerate the app icon
```sh
python assets/gen_icon.py            # -> assets/icon.png (window) + assets/icon.ico (exe/installer)
```

## Linux

The stack is pure Rust with no platform-specific code outside `build.rs` (Windows-only,
cfg-gated). `rfd` uses the XDG desktop portal (no GTK) and `ureq` uses rustls (no OpenSSL),
so the build deps are minimal:

```sh
sudo apt-get install -y libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libxkbcommon-dev
cargo build --release            # -> target/release/GenomeForge
```

Desktop integration files (used by the tarball and Flatpak) live in `packaging/linux/`:
the `.desktop` entry and the AppStream `.metainfo.xml`, both named with the app ID
`io.github.abowlofeleven.GenomeForge`.

### Flatpak
Manifest: `packaging/flatpak/io.github.abowlofeleven.GenomeForge.yml`. It builds the app
**inside** the freedesktop SDK (24.08) with **network access**: the build installs the
pinned Rust toolchain (1.96.0) via `rustup` and lets cargo fetch crates from crates.io.
(The SDK's bundled `rust-stable` extension is 1.89, too old for egui 0.34's 1.92 MSRV, so
that route is not used.) No vendoring or `cargo-sources.json` is needed:

```sh
flatpak install -y flathub org.freedesktop.Platform//24.08 org.freedesktop.Sdk//24.08
flatpak-builder --user --install --force-clean build-dir \
  packaging/flatpak/io.github.abowlofeleven.GenomeForge.yml
flatpak run io.github.abowlofeleven.GenomeForge
```
File access is portal-mediated, so the sandbox needs no broad `--filesystem` permission.

> The network build is simple and fine for a self-hosted bundle, but **not Flathub-compliant**
> (Flathub requires fully offline builds). Submitting there later would mean switching back to
> vendored sources with a newer toolchain extension. The CI release job marks the Flatpak step
> `continue-on-error`, so a Flatpak hiccup never blocks a release.

## macOS

```sh
# Universal (Intel + Apple Silicon):
rustup target add x86_64-apple-darwin aarch64-apple-darwin
cargo build --release --target x86_64-apple-darwin -p gx-app
cargo build --release --target aarch64-apple-darwin -p gx-app
lipo -create -output GenomeForge \
  target/x86_64-apple-darwin/release/GenomeForge target/aarch64-apple-darwin/release/GenomeForge
```
`packaging/macos/Info.plist` is the bundle plist; the CI release job assembles
`GenomeForge.app` (with an `.icns` built from `assets/icon.png`), ad-hoc signs it
(`codesign -s -`, so the arm64 slice launches), and wraps it in a `.dmg`. Builds are
**unsigned/unnotarized**, so Gatekeeper will warn; right-click ▸ Open (or
`xattr -dr com.apple.quarantine GenomeForge.app`) to run.

## Continuous integration & releases

- **`.github/workflows/ci.yml`**: builds, tests, and runs `clippy -D warnings` on
  Windows, Linux, and macOS for every push to `main` and every PR.
- **`.github/workflows/release.yml`**: on a `v*` tag (or the manual "Run workflow"
  button), re-verifies on each OS and, only if that passes, builds and publishes a single
  GitHub Release with: Windows `.exe` + `.msi`, Linux `.tar.gz` + Flatpak bundle, and a
  universal macOS `.dmg`. Cut a release with:

  ```sh
  git tag v0.1.0 && git push origin v0.1.0
  ```

  Bumping the version means updating **all four** version sources before tagging:
  `Cargo.toml` (`[workspace.package].version`), `installer/genomeforge.wxs` (`Version`),
  `packaging/macos/Info.plist`, and `packaging/linux/*.metainfo.xml`. The release job
  fails fast if the tag and `Cargo.toml` disagree.
