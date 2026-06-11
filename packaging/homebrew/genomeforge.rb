# Homebrew cask for GenomeForge.
#
# Until the macOS build is notarized, distribute this via a personal tap rather
# than homebrew-cask core (which requires signed/notarized apps):
#
#   brew tap ABowlOfEleven/tap https://github.com/ABowlOfEleven/homebrew-tap
#   brew install --cask genomeforge
#
# Bump `version` + `sha256` on each release (sha256 of the .dmg).
cask "genomeforge" do
  version "0.1.2"
  sha256 "7395a4a30c6c2d1c8bf5fe75eaab22ccb96995832da2e3101c97227f2a519836"

  url "https://github.com/ABowlOfEleven/genomeforge/releases/download/v#{version}/GenomeForge-#{version}-macos-universal.dmg",
      verified: "github.com/ABowlOfEleven/genomeforge/"
  name "GenomeForge"
  desc "Local genome browser, variant explorer, plasmid designer, and CRISPR studio"
  homepage "https://github.com/ABowlOfEleven/genomeforge"

  app "GenomeForge.app"

  # The build is currently unsigned/unnotarized; clear the quarantine flag so it
  # launches without the Gatekeeper "damaged" prompt.
  postflight do
    system_command "/usr/bin/xattr",
                   args: ["-dr", "com.apple.quarantine", "#{appdir}/GenomeForge.app"]
  end

  zap trash: [
    "~/Library/Application Support/GenomeForge",
    "~/Library/Caches/GenomeForge",
  ]
end
