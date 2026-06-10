# winget manifest

Draft manifest for submitting GenomeForge to the
[Windows Package Manager Community Repository](https://github.com/microsoft/winget-pkgs).

The three YAML files here go into `manifests/a/ABowlOfEleven/GenomeForge/<version>/` in a
fork of `microsoft/winget-pkgs`, then opened as a PR. The easiest path is
[`wingetcreate`](https://github.com/microsoft/winget-create):

```pwsh
winget install wingetcreate
wingetcreate update ABowlOfEleven.GenomeForge `
  --version 0.1.1 `
  --urls https://github.com/ABowlOfEleven/genomeforge/releases/download/v0.1.1/GenomeForge-0.1.1-x64.msi `
  --submit
```

On a new release, bump `PackageVersion` in all three files and the installer URL +
`InstallerSha256` (uppercase, no `sha256:` prefix) in the installer manifest.
