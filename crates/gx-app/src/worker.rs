//! Background worker thread that owns the [`AnnotationService`] (and its redb
//! cache) and serves the UI over channels. Keeps all network + disk I/O off the
//! paint thread; calls `ctx.request_repaint()` after each response so the UI
//! refreshes promptly.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::mpsc::{Receiver, Sender, channel};
use std::thread;

use eframe::egui;

use gx_core::{Assembly, Feature, GenomicRange};
use gx_annotate::{
    AnnotationService, Article, Cache, CacheStats, GeneLocation, UpdateInfo, VariantAnnotation,
};

pub enum Request {
    Import(PathBuf),
    Annotate(Vec<String>),
    Features(GenomicRange),
    Sequence(GenomicRange),
    LookupGene(String),
    FetchPgs(String),
    /// Fetch related PubMed articles for a variant (by rsID).
    Literature(String),
    /// Lift a coordinate from one assembly to another.
    Liftover {
        from: Assembly,
        to: Assembly,
        contig: String,
        pos: u64,
    },
    /// Check GitHub for a newer release of the app.
    CheckUpdate(String),
    /// Report how many entries are in the local annotation cache.
    CacheStats,
    /// Empty the local annotation cache.
    ClearCache,
    SetOnline(bool),
    SetAssembly(Assembly),
    Shutdown,
}

#[allow(dead_code)] // `Notice` and `Features.range` are wired for future use.
pub enum Response {
    Imported(Result<gx_io::Imported, String>),
    Annotations(HashMap<String, VariantAnnotation>),
    Features {
        range: GenomicRange,
        features: Vec<Feature>,
    },
    Sequence {
        range: GenomicRange,
        seq: String,
    },
    Gene {
        symbol: String,
        location: Option<GeneLocation>,
    },
    Pgs(Result<gx_pgs::ScoreFile, String>),
    /// PubMed articles for a variant (keyed by the rsID that was searched).
    Literature {
        rsid: String,
        articles: Vec<Article>,
    },
    /// Liftover result: the lifted `(contig, pos)`, or `None` if unmapped.
    Liftover {
        to: Assembly,
        mapped: Option<(String, u64)>,
    },
    /// Latest GitHub release, if the check succeeded (`None` = up to date / no releases).
    Update(Option<UpdateInfo>),
    /// Current (or post-clear) local cache entry counts.
    CacheStats(CacheStats),
    Notice(String),
    Failed(String),
    /// A counted request that produced no data (keeps the in-flight tally exact).
    Idle,
}

pub struct Worker {
    tx: Sender<Request>,
    pub rx: Receiver<Response>,
}

impl Worker {
    pub fn spawn(
        ctx: egui::Context,
        cache_path: PathBuf,
        online: bool,
        api_key: Option<String>,
    ) -> Self {
        let (req_tx, req_rx) = channel::<Request>();
        let (res_tx, res_rx) = channel::<Response>();

        thread::spawn(move || {
            let cache = match Cache::open(&cache_path) {
                Ok(c) => c,
                Err(e) => {
                    let _ = res_tx.send(Response::Failed(format!("cache init failed: {e}")));
                    return;
                }
            };
            let mut service = AnnotationService::new(cache, online, Assembly::Grch38, api_key);

            while let Ok(req) = req_rx.recv() {
                let response = handle(&mut service, req);
                match response {
                    HandleOutcome::Reply(r) => {
                        let _ = res_tx.send(r);
                        ctx.request_repaint();
                    }
                    HandleOutcome::Silent => {}
                    HandleOutcome::Stop => break,
                }
            }
        });

        Self {
            tx: req_tx,
            rx: res_rx,
        }
    }

    pub fn send(&self, req: Request) {
        let _ = self.tx.send(req);
    }
}

impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.tx.send(Request::Shutdown);
    }
}

enum HandleOutcome {
    Reply(Response),
    Silent,
    Stop,
}

fn handle(service: &mut AnnotationService, req: Request) -> HandleOutcome {
    match req {
        Request::Shutdown => HandleOutcome::Stop,
        Request::SetOnline(b) => {
            service.set_online(b);
            HandleOutcome::Silent
        }
        Request::SetAssembly(a) => {
            service.set_assembly(a);
            HandleOutcome::Silent
        }
        Request::Import(path) => {
            log::info!("importing {}", path.display());
            let result = gx_io::import_path(&path).map_err(|e| e.to_string());
            match &result {
                Ok(gx_io::Imported::Variants(v)) => {
                    log::info!("imported {} variants ({})", v.store.len(), v.format_label)
                }
                Ok(gx_io::Imported::Sequences(s)) => log::info!("imported {} sequence(s)", s.len()),
                Err(e) => log::warn!("import failed: {e}"),
            }
            HandleOutcome::Reply(Response::Imported(result))
        }
        Request::Annotate(rsids) => match service.annotate(&rsids) {
            Ok(map) => {
                log::info!("annotated {} of {} requested rsIDs", map.len(), rsids.len());
                HandleOutcome::Reply(Response::Annotations(map))
            }
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("annotation: {e}"))),
        },
        Request::Features(range) => match service.features(&range) {
            Ok(features) => {
                log::info!("{} features in {}", features.len(), range.ensembl_region());
                HandleOutcome::Reply(Response::Features { range, features })
            }
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("features: {e}"))),
        },
        Request::Sequence(range) => match service.sequence(&range) {
            Ok(Some(seq)) => {
                log::info!("sequence {} ({} bp)", range.ensembl_region(), seq.len());
                HandleOutcome::Reply(Response::Sequence { range, seq })
            }
            // Counted request: reply (even with no data) so the in-flight tally stays exact.
            Ok(None) => HandleOutcome::Reply(Response::Idle),
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("sequence: {e}"))),
        },
        Request::LookupGene(symbol) => match service.gene_location(&symbol) {
            Ok(location) => {
                log::info!("gene {symbol}: {}", if location.is_some() { "found" } else { "not found" });
                HandleOutcome::Reply(Response::Gene { symbol, location })
            }
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("gene lookup: {e}"))),
        },
        Request::FetchPgs(id) => match service.fetch_pgs_text(&id) {
            Ok(text) => {
                log::info!("fetched PGS {id} ({} bytes)", text.len());
                HandleOutcome::Reply(Response::Pgs(gx_pgs::parse(&text)))
            }
            Err(e) => HandleOutcome::Reply(Response::Pgs(Err(format!("download: {e}")))),
        },
        Request::Literature(rsid) => match service.literature(&rsid, 6) {
            Ok(articles) => {
                log::info!("PubMed {rsid}: {} article(s)", articles.len());
                HandleOutcome::Reply(Response::Literature { rsid, articles })
            }
            Err(e) => {
                log::warn!("PubMed {rsid} failed: {e}");
                // Reply with an empty set so the UI clears "searching…" and the
                // in-flight tally stays exact.
                HandleOutcome::Reply(Response::Literature { rsid, articles: Vec::new() })
            }
        },
        Request::Liftover { from, to, contig, pos } => match service.liftover(from, to, &contig, pos) {
            Ok(mapped) => {
                log::info!("liftover {contig}:{pos} {} -> {}: {mapped:?}", from.label(), to.label());
                HandleOutcome::Reply(Response::Liftover { to, mapped })
            }
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("liftover: {e}"))),
        },
        Request::CheckUpdate(repo) => {
            if !service.is_online() {
                return HandleOutcome::Reply(Response::Idle);
            }
            match gx_annotate::latest_release(&repo) {
                Ok(info) => HandleOutcome::Reply(Response::Update(info)),
                Err(e) => {
                    log::warn!("update check failed: {e}");
                    HandleOutcome::Reply(Response::Update(None))
                }
            }
        }
        Request::CacheStats => match service.cache_stats() {
            Ok(s) => HandleOutcome::Reply(Response::CacheStats(s)),
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("cache stats: {e}"))),
        },
        Request::ClearCache => match service.clear_cache() {
            Ok(()) => {
                log::info!("annotation cache cleared");
                let stats = service.cache_stats().unwrap_or_default();
                HandleOutcome::Reply(Response::CacheStats(stats))
            }
            Err(e) => HandleOutcome::Reply(Response::Failed(format!("clear cache: {e}"))),
        },
    }
}
