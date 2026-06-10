// Hide the console window in release builds (keep it in debug for logs).
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod browser;
mod crispr;
mod document;
mod health;
mod phenotype;
mod plasmid;
mod settings;
mod theme;
mod tutorial;
mod ux;
mod view;
mod worker;

use eframe::egui;

use app::GenomeForgeApp;

fn main() -> eframe::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let mut viewport = egui::ViewportBuilder::default()
        .with_inner_size([1280.0, 800.0])
        .with_min_inner_size([860.0, 540.0])
        .with_title("GenomeForge");
    if let Ok(icon) = eframe::icon_data::from_png_bytes(include_bytes!("../../../assets/icon.png")) {
        viewport = viewport.with_icon(std::sync::Arc::new(icon));
    }
    let native_options = eframe::NativeOptions {
        viewport,
        ..Default::default()
    };

    eframe::run_native(
        "GenomeForge",
        native_options,
        Box::new(|cc| Ok(Box::new(GenomeForgeApp::new(cc)))),
    )
}
