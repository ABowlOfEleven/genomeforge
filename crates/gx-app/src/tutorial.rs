//! In-app tutorials. Each workspace has a short, navigable walkthrough that can
//! be opened from the Help button and auto-shows the first time a workspace is
//! entered.

use eframe::egui;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Section {
    Genome,
    Plasmid,
    Crispr,
    Phenotype,
}

pub struct Step {
    pub title: &'static str,
    pub body: &'static str,
}

pub fn title(section: Section) -> &'static str {
    match section {
        Section::Genome => "Genome browser",
        Section::Plasmid => "Plasmid designer",
        Section::Crispr => "CRISPR studio",
        Section::Phenotype => "Phenotype & risk",
    }
}

pub fn steps(section: Section) -> &'static [Step] {
    match section {
        Section::Genome => GENOME,
        Section::Plasmid => PLASMID,
        Section::Crispr => CRISPR,
        Section::Phenotype => PHENOTYPE,
    }
}

const GENOME: &[Step] = &[
    Step {
        title: "Welcome to the genome browser",
        body: "This is where you explore a genome — yours or anyone's. Open a file with \
               File ▸ Open: a 23andMe/AncestryDNA download, a VCF from sequencing, or a \
               FASTA/GenBank sequence. Your raw data never leaves your computer.",
    },
    Step {
        title: "Reading the tracks",
        body: "From top to bottom: a ruler (genomic position), an ideogram (the whole \
               chromosome, with your current window highlighted), gene models (boxes are \
               exons, lines are introns), and your variants drawn as lollipops. Drag to pan, \
               scroll to zoom. Zoom all the way in to see the DNA letters.",
    },
    Step {
        title: "Variants and colour",
        body: "Each lollipop is a position where your DNA was read. Colour shows its ClinVar \
               clinical significance — red = pathogenic, green = benign, yellow = uncertain. \
               Click one to see its full annotation on the right.",
    },
    Step {
        title: "Search and the table",
        body: "Use the search box for a gene name (BRCA1), an rsID (rs6025), or a position \
               (17:43,044,295). The sidebar's 'Common health SNPs' list jumps straight to \
               well-known variants (APOE, Factor V Leiden, MTHFR…). The variant table at the \
               bottom is filterable — try 'Notable only' to see clinically meaningful variants.",
    },
    Step {
        title: "A note on the science",
        body: "Annotations come from public databases (ClinVar, gnomAD, dbSNP) via \
               MyVariant.info, cached locally. This is a research & learning tool — not \
               medical advice. Talk to a clinician or genetic counsellor about real results.",
    },
];

const PLASMID: &[Step] = &[
    Step {
        title: "Welcome to the plasmid designer",
        body: "A plasmid is a small, usually circular piece of DNA used to carry genes into \
               cells. Open a FASTA or GenBank file and you land here. The big circle is a map \
               of the whole molecule.",
    },
    Step {
        title: "The circular & linear maps",
        body: "Coloured arcs on the circle are features (genes, promoters, origins). Ticks \
               outside the ring are restriction-enzyme cut sites. Below is a linear view of \
               the same sequence. Drag the circle to rotate it; pan/zoom the linear view.",
    },
    Step {
        title: "Restriction sites & unique cutters",
        body: "Restriction enzymes cut DNA at specific sequences. 'Unique cutters' (in the \
               sidebar) cut the plasmid exactly once — those are the handy ones for cloning, \
               because they give you a single, predictable place to insert DNA.",
    },
    Step {
        title: "ORFs, translation & primers",
        body: "ORFs are stretches that could encode protein (ATG…stop). Click a feature or \
               ORF to translate it and read the protein. Selecting a feature also shows its \
               melting temperature and GC content — useful when designing primers.",
    },
];

const CRISPR: &[Step] = &[
    Step {
        title: "Welcome to CRISPR studio",
        body: "CRISPR-Cas9 is molecular scissors: a short 'guide RNA' leads the Cas9 protein \
               to a matching 20-letter stretch of DNA and cuts it. Here you design guides for \
               a target sequence and explore the edits they'd make.",
    },
    Step {
        title: "Choosing a target",
        body: "In the sidebar, click 'Use loaded sequence' to design against the plasmid/gene \
               you already opened, or paste any DNA and press Analyze. The studio finds every \
               valid guide site (those followed by an 'NGG' PAM) on both strands.",
    },
    Step {
        title: "Reading the guide table",
        body: "Each row is a candidate guide. The on-target estimate suggests how efficiently \
               it should cut; GC content and a poly-T warning flag guides that may work poorly. \
               Click a guide to inspect it.",
    },
    Step {
        title: "Off-targets & specificity",
        body: "Cas9 can also cut places that *almost* match the guide. The studio scores \
               specificity — mismatches near the PAM (the 'seed') matter most. It searches the \
               target by default; load a reference FASTA (a chromosome or whole genome) under \
               'Off-target reference' to check genome-wide.",
    },
    Step {
        title: "Edits & AAV delivery",
        body: "Simulate the repair outcome: a small NHEJ deletion (often disabling a gene) or a \
               precise HDR replacement with a donor. The AAV cargo planner checks whether your \
               construct fits inside an AAV virus (~4.7 kb), a common delivery vehicle.",
    },
    Step {
        title: "Use responsibly",
        body: "Scores here are transparent decision-support estimates, not a clinical pipeline. \
               Genome editing of humans is tightly regulated — this tool is for learning, \
               research, and design exploration.",
    },
];

const PHENOTYPE: &[Step] = &[
    Step {
        title: "Phenotype & risk",
        body: "Most traits — height, common disease risk — are shaped by thousands of genetic \
               variants, each with a tiny effect. A polygenic score (PRS) adds them up to estimate \
               your genetic tendency. It's a probability, never a diagnosis.",
    },
    Step {
        title: "Fetching a score",
        body: "Pick an example score or enter any PGS Catalog ID (like PGS000001), or load a \
               scoring file from disk. GenomeForge downloads it, matches its variants to your \
               genome, and computes your score — handling strand flips automatically.",
    },
    Step {
        title: "Reading the result",
        body: "'Coverage' is how many of the score's variants your data actually contains — \
               consumer chips cover less than whole-genome sequencing, so watch for low coverage. \
               When the file includes population frequencies, you also get an estimated percentile.",
    },
    Step {
        title: "Big caveats",
        body: "Percentiles use a normal approximation and assume the reference population is like \
               you — most scores were trained on European-ancestry data and transfer poorly to \
               others. Environment and lifestyle usually matter as much as genetics.",
    },
    Step {
        title: "Mendelian findings",
        body: "On the right, single variants ClinVar flags as clinically significant. Unlike a \
               PRS, one of these can matter on its own. None of this is medical advice — take real \
               questions to a clinician or genetic counsellor.",
    },
];

/// Draw the tutorial window. Sets `*open = false` when dismissed.
pub fn show(ctx: &egui::Context, section: Section, step: &mut usize, open: &mut bool) {
    let steps = steps(section);
    let n = steps.len().max(1);
    *step = (*step).min(n - 1);

    let mut keep_open = true;
    egui::Window::new(format!("Tutorial — {}", title(section)))
        .open(&mut keep_open)
        .collapsible(false)
        .resizable(false)
        .default_width(440.0)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            let s = &steps[*step];
            ui.heading(s.title);
            ui.add_space(6.0);
            ui.label(egui::RichText::new(s.body).size(13.0));
            ui.add_space(12.0);
            ui.separator();
            ui.horizontal(|ui| {
                if ui.add_enabled(*step > 0, egui::Button::new("◀ Back")).clicked() {
                    *step -= 1;
                }
                ui.label(format!("{} / {}", *step + 1, n));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if *step + 1 >= n {
                        if ui.button("Done").clicked() {
                            *open = false;
                        }
                    } else if ui.button("Next ▶").clicked() {
                        *step += 1;
                    }
                });
            });
        });

    if !keep_open {
        *open = false;
    }
}
