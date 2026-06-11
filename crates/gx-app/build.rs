//! Build script: embed the app icon + version metadata into the Windows .exe.
//! Fails soft — if no resource compiler is available, the build still succeeds
//! (just without the embedded icon).

fn main() {
    println!("cargo:rerun-if-changed=../../assets/icon.ico");

    #[cfg(windows)]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("../../assets/icon.ico");
        res.set("ProductName", "GenomeForge");
        res.set(
            "FileDescription",
            "GenomeForge: genome browser & variant explorer",
        );
        res.set("CompanyName", "GenomeForge");
        res.set("LegalCopyright", "GenomeForge contributors (MIT)");
        if let Err(e) = res.compile() {
            println!("cargo:warning=icon/metadata embed skipped: {e}");
        }
    }
}
