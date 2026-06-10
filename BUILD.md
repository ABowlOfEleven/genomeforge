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
> Note: WiX v6/v7 require accepting a paid "Open Source Maintenance Fee" EULA — stay on v5.

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

## Cross-compiling to Linux (future)
The stack is pure Rust with no platform-specific code outside `build.rs` (Windows-only,
cfg-gated), so a Linux build is `cargo build --release` on Linux (or a cross toolchain).
A `.deb`/AppImage can be added later with `cargo-deb` / `cargo-bundle`.
