// On Windows the app is a GUI program: no console window opens. Output still reaches a redirected file.
#![cfg_attr(windows, windows_subsystem = "windows")]
//! FoldShift: local server (127.0.0.1 only) with a one-page UI. Exits when the page is closed.

use psplot_core::compare::{
    compare, distance_table_csv, marked_rows, plain_summary, substituted_residues, Comparison,
    Options,
};
use psplot_core::structure::{load_any, transform_cif, transform_pdb};
use psplot_core::transcript::{from_sequence, get_transcript};
use psplot_core::variants::apply_variants;
use psplot_render::{
    encode_jpeg, figure, profile_figure, view_data, Kind, Style, PREVIEW_DPI, PRINT_DPI,
};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::io::Cursor;
use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use tiny_http::{Header, Method, Request, Response, Server};

const INDEX: &str = include_str!("../index.html");
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A finished comparison kept in memory so the rotating view is quick.
type Cached = Arc<(Comparison, Params)>;

struct State {
    runs: PathBuf,
    server: Option<ServerCfg>,
    cache: PathBuf,
    jobs: Mutex<HashMap<String, Value>>,
    busy: AtomicUsize,
    last_ping: AtomicI64, // ms since start
    bye_at: AtomicI64,    // ms since start, 0 = none
    seen: AtomicBool,
    start: Instant,
    counter: AtomicUsize,
    cache_cmp: Mutex<Vec<(String, Cached)>>,
}

/// Settings of the shared server mode (`--server` or `PSPLOT_SERVER=1`), for use in a container.
struct ServerCfg {
    password: Option<String>,
    max_upload: usize,
    max_jobs: usize,
    keep_days: u64,
}

thread_local! {
    /// Browser session of the request being handled (server mode only).
    static SID: std::cell::RefCell<String> = const { std::cell::RefCell::new(String::new()) };
}

impl State {
    /// Folder with the runs of the current request: the shared folder on a desktop, one folder per browser session on a server.
    fn dir(&self) -> PathBuf {
        if self.server.is_some() {
            SID.with(|s| self.runs.join(s.borrow().as_str()))
        } else {
            self.runs.clone()
        }
    }
}

fn random_hex(n: usize) -> String {
    use std::io::Read as _;
    let mut b = vec![0u8; n];
    let ok = std::fs::File::open("/dev/urandom")
        .and_then(|mut f| f.read_exact(&mut b))
        .is_ok();
    if !ok {
        let t = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        for (i, x) in b.iter_mut().enumerate() {
            *x = ((t >> ((i % 16) * 8)) as u8)
                ^ (i as u8).wrapping_mul(131)
                ^ std::process::id() as u8;
        }
    }
    b.iter().map(|x| format!("{x:02x}")).collect()
}

fn b64_decode(t: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let (mut acc, mut bits) = (0u32, 0);
    for c in t.bytes() {
        let v = match c {
            b'A'..=b'Z' => c - b'A',
            b'a'..=b'z' => c - b'a' + 26,
            b'0'..=b'9' => c - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            _ => continue,
        };
        acc = (acc << 6) | v as u32;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((acc >> bits) as u8);
            acc &= (1 << bits) - 1;
        }
    }
    out
}

fn same_secret(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).fold(0u8, |d, (x, y)| d | (x ^ y)) == 0
}

fn now_ms(st: &State) -> i64 {
    st.start.elapsed().as_millis() as i64
}

fn data_home() -> PathBuf {
    if let Ok(p) = std::env::var("PSPLOT_HOME") {
        return PathBuf::from(p);
    }
    let home = |v: &str| std::env::var(v).map(PathBuf::from).ok();
    if cfg!(windows) {
        home("LOCALAPPDATA")
            .or_else(|| home("USERPROFILE"))
            .unwrap_or_else(|| PathBuf::from("."))
            .join("FoldShift")
    } else if cfg!(target_os = "macos") {
        home("HOME")
            .unwrap_or_else(|| PathBuf::from("."))
            .join("Library/Application Support/FoldShift")
    } else {
        home("XDG_DATA_HOME")
            .or_else(|| home("HOME").map(|h| h.join(".local/share")))
            .unwrap_or_else(|| PathBuf::from("."))
            .join("FoldShift")
    }
}

fn json_resp(code: u16, v: Value) -> Response<Cursor<Vec<u8>>> {
    Response::from_string(v.to_string())
        .with_status_code(code)
        .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
}

fn err(code: u16, m: &str) -> Response<Cursor<Vec<u8>>> {
    json_resp(code, json!({ "error": m }))
}

// ---------------------------------------------------------------- multipart
struct Part {
    name: String,
    data: Vec<u8>,
}

fn find(h: &[u8], n: &[u8], from: usize) -> Option<usize> {
    h.get(from..)?
        .windows(n.len())
        .position(|w| w == n)
        .map(|p| p + from)
}

fn multipart(body: &[u8], ctype: &str) -> Vec<Part> {
    let Some(b) = ctype.split("boundary=").nth(1) else {
        return vec![];
    };
    let delim = format!("--{}", b.trim_matches('"')).into_bytes();
    let mut parts = Vec::new();
    let mut pos = find(body, &delim, 0);
    while let Some(p) = pos {
        let start = p + delim.len();
        if body.get(start..start + 2) == Some(b"--") {
            break;
        }
        let Some(hdr_end) = find(body, b"\r\n\r\n", start) else {
            break;
        };
        let headers = String::from_utf8_lossy(&body[start..hdr_end]).to_string();
        let next = find(body, &delim, hdr_end + 4).unwrap_or(body.len());
        let data_end = next.saturating_sub(2).max(hdr_end + 4);
        let name = headers
            .split("name=\"")
            .nth(1)
            .and_then(|s| s.split('"').next())
            .unwrap_or("")
            .to_string();
        parts.push(Part {
            name,
            data: body[hdr_end + 4..data_end].to_vec(),
        });
        pos = Some(next).filter(|n| *n < body.len());
    }
    parts
}

// ---------------------------------------------------------------- handlers
fn sequence(st: &State, body: &[u8]) -> Response<Cursor<Vec<u8>>> {
    let Ok(d) = serde_json::from_slice::<Value>(body) else {
        return err(400, "Bad request.");
    };
    let s = |k: &str| d[k].as_str().unwrap_or("").to_string();
    let tx = if !s("sequence").trim().is_empty() {
        from_sequence(&s("sequence"))
    } else {
        get_transcript(&s("transcript"), &st.cache)
    };
    let tx = match tx {
        Ok(t) => t,
        Err(e) => return err(400, &e.0),
    };
    let vs: Vec<String> = s("variants")
        .split(['\n', ',', ';'])
        .map(|x| x.trim().to_string())
        .filter(|x| !x.is_empty())
        .collect();
    let res = match apply_variants(&tx, &vs) {
        Ok(r) => r,
        Err(e) => return err(400, &e.0),
    };
    let gene = if tx.gene.is_empty() {
        tx.accession.clone()
    } else {
        tx.gene.clone()
    };
    let ref_name = format!("{gene} reference {}", tx.accession);
    let case_name = format!(
        "{gene} case {}",
        if res.variants.is_empty() {
            "no variant".to_string()
        } else {
            res.variants.join(" ")
        }
    );
    json_resp(
        200,
        json!({
            "gene": gene, "accession": tx.accession, "definition": tx.definition,
            "reference": tx.protein, "case": res.protein, "changed": res.changed, "notes": res.notes,
            "variants": res.variants, "marks": substituted_residues(&tx.protein, &res.protein).iter().take(30).map(|x| x.to_string()).collect::<Vec<_>>().join(", "), "ref_name": ref_name, "case_name": case_name
        }),
    )
}

struct Params {
    opt: Options,
    style: Style,
    /// Marked residues are also the focus when `zoom_on_marks` is set.
    raw: HashMap<String, String>,
}

/// Residue numbers typed by the user (marks, zoom, alignment region) are numbers of the reference file, such as 508 for
/// Phe508. Turn them into positions 1, 2, 3, ... of the loaded reference. Numbers that are not in the file are dropped.
/// Files whose numbers repeat (several chains) keep the position meaning.
fn apply_numbering(p: &mut Params, r: &psplot_core::structure::AfResult) {
    let nums = psplot_core::compare::numbering(&r.models[r.best()].resnums);
    if nums.iter().enumerate().all(|(i, &n)| n == i as i32 + 1) {
        return;
    }
    let pos: HashMap<i64, usize> = nums
        .iter()
        .enumerate()
        .map(|(i, &n)| (n as i64, i + 1))
        .collect();
    let map = |v: &[usize]| -> Vec<usize> {
        v.iter()
            .filter_map(|&n| pos.get(&(n as i64)).copied())
            .collect()
    };
    p.style.marks = map(&p.style.marks);
    p.style.focus = map(&p.style.focus);
    if !p.opt.fit_region.trim().is_empty() {
        p.opt.fit_region = map(&residue_list(&p.opt.fit_region))
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(",");
        if p.opt.fit_region.is_empty() {
            p.opt.fit_region = "0".into(); // none of the typed residues exists: an empty selection, not "automatic"
        }
    }
}

fn residue_list(text: &str) -> Vec<usize> {
    // "12, 20-25" -> [12, 20, 21, ..., 25]
    let mut out = Vec::new();
    for part in text.split([',', ';', ' ', '\n']).filter(|p| !p.is_empty()) {
        let nums: Vec<usize> = part
            .split('-')
            .filter_map(|x| x.trim().parse().ok())
            .collect();
        match nums.as_slice() {
            [a] => out.push(*a),
            [a, b] if a <= b && b - a < 5000 => out.extend(*a..=*b),
            _ => {}
        }
    }
    out
}

fn parse_params(f: &HashMap<String, String>) -> Result<Params, String> {
    let num = |k: &str, d: f64| -> Result<f64, String> {
        match f.get(k).map(|s| s.trim()).filter(|s| !s.is_empty()) {
            None => Ok(d),
            Some(s) => s
                .parse::<f64>()
                .map_err(|_| "An option is not a number.".to_string()),
        }
    };
    let text = |k: &str| f.get(k).cloned().unwrap_or_default();
    let style = Style {
        cap: num("cap", 10.0)?.max(1.0),
        view: num("view", 0.0)? as u8,
        marks: residue_list(&text("mark")),
        focus: residue_list(&text("focus")),
        plddt_panel: f.get("plddt").map(|v| v.trim()) != Some("0"),
        yaw: num("yaw", 0.0)? as f32,
        pitch: num("pitch", 0.0)? as f32,
        roll: num("roll", 0.0)? as f32,
        zoom: num("zoom", 1.0)? as f32,
        steady: false,
    };
    Ok(Params {
        opt: Options {
            fit_region: text("fit_region"),
            min_plddt: num("min_plddt", 70.0)?,
            cutoff: num("cutoff", 2.0)?,
            robustness: true,
            blosum: f.get("pairing").map(|v| v.trim()) != Some("identity"),
        },
        style,
        raw: f.clone(),
    })
}

fn write(dir: &Path, name: &str, data: &[u8]) -> Result<(), String> {
    std::fs::write(dir.join(name), data).map_err(|e| format!("Could not write {name}: {e}"))
}

/// Load the saved upload of one side ("ref" or "case"): `<side>.zip` holds the structure file or zip
/// (whatever its type), `<side>.name` its original name, `<side>.pae` an optional PAE file.
fn load_side(
    dir: &Path,
    side: &str,
    label: &str,
) -> Result<psplot_core::structure::AfResult, String> {
    let bytes = std::fs::read(dir.join(format!("{side}.zip")))
        .map_err(|e| format!("Could not open the {label} file: {e}"))?;
    let name = std::fs::read_to_string(dir.join(format!("{side}.name"))).unwrap_or_default();
    let pae = std::fs::read(dir.join(format!("{side}.pae"))).ok();
    load_any(&bytes, &name, pae.as_deref()).map_err(|e| format!("{label} file: {e}"))
}

/// Read one side and keep the chosen chains.
fn load_one(
    dir: &Path,
    side: &str,
    label: &str,
    chains: &str,
) -> Result<psplot_core::structure::AfResult, String> {
    let mut a = load_side(dir, side, label)?;
    a.keep_chains(chains)
        .map_err(|e| format!("{label} file: {e}"))?;
    Ok(a)
}

/// Read the reference and the case at the same time.
fn load_pair(
    dir: &Path,
    chains: &str,
) -> Result<
    (
        psplot_core::structure::AfResult,
        psplot_core::structure::AfResult,
    ),
    String,
> {
    std::thread::scope(|sc| {
        let r = sc.spawn(|| load_one(dir, "ref", "Reference", chains));
        let c = load_one(dir, "case", "Case", chains);
        let r = r
            .join()
            .unwrap_or_else(|_| Err("Internal error while reading the file.".into()));
        Ok((r?, c?))
    })
}

const SIDE_FILES: [&str; 6] = [
    "ref.zip",
    "ref.name",
    "ref.pae",
    "case.zip",
    "case.name",
    "case.pae",
];

fn run_job(st: &State, id: &str, dir: &Path, mut p: Params) -> Result<Value, String> {
    let (r, c) = load_pair(dir, p.raw.get("chains").map_or("", |s| s.as_str()))?;
    apply_numbering(&mut p, &r);
    let cmp = compare(&r, &c, &p.opt)?;
    outputs(dir, &cmp, &p)?;
    let mut raw = p.raw.clone();
    raw.remove("from_job");
    raw.remove("ref_from");
    write(
        dir,
        "params.json",
        serde_json::to_string(&raw).unwrap().as_bytes(),
    )?;
    let plain = plain_summary(&cmp);
    write(dir, "summary.txt", plain.join("\n").as_bytes())?;
    let summary = json!({
        "stats": cmp.stats, "blocks": cmp.blocks, "robustness": cmp.robustness,
        "missing_in_case": cmp.missing_in_case.iter().map(|&(a, b)| (cmp.ref_num[a - 1], cmp.ref_num[b - 1])).collect::<Vec<_>>(), "extra_in_case": cmp.extra_in_case.iter().map(|&(a, b)| (cmp.case_num[a - 1], cmp.case_num[b - 1])).collect::<Vec<_>>(), "plain": plain,
        "ss_changes": cmp.ss_changes, "marked": marked_rows(&cmp, &p.style.marks),
        "reference": cmp.ref_name, "case": cmp.case_name, "version": VERSION,
        "has_plddt": p.style.plddt_panel && cmp.ref_model.has_conf && cmp.case_model.has_conf,
        "has_pae": cmp.ref_model.pae.is_some() && cmp.case_model.pae.is_some(),
        "sources": [cmp.ref_source, cmp.case_source]
    });
    write(
        dir,
        "summary.json",
        serde_json::to_string_pretty(&summary).unwrap().as_bytes(),
    )?;
    remember(st, id, cmp, p);
    Ok(summary)
}

fn remember(st: &State, id: &str, cmp: Comparison, p: Params) {
    let mut c = st.cache_cmp.lock().unwrap();
    c.retain(|(k, _)| k != id);
    c.push((id.to_string(), Arc::new((cmp, p))));
    while c.len() > 4 {
        c.remove(0);
    }
}

/// The comparison of a run: from memory, or rebuilt from the stored zips and parameters.
fn comparison_of(st: &State, id: &str) -> Result<Cached, String> {
    if let Some((_, v)) = st.cache_cmp.lock().unwrap().iter().find(|(k, _)| k == id) {
        return Ok(Arc::clone(v));
    }
    if !valid_id(id) {
        return Err("Unknown run.".into());
    }
    let dir = st.dir().join(id);
    let raw: HashMap<String, String> = std::fs::read_to_string(dir.join("params.json"))
        .ok()
        .and_then(|t| serde_json::from_str(&t).ok())
        .ok_or("This run is no longer available.")?;
    let mut p = parse_params(&raw)?;
    let (r, c) = load_pair(&dir, p.raw.get("chains").map_or("", |s| s.as_str()))?;
    apply_numbering(&mut p, &r);
    let cmp = compare(&r, &c, &p.opt)?;
    let v = Arc::new((cmp, p));
    let mut cc = st.cache_cmp.lock().unwrap();
    cc.push((id.to_string(), Arc::clone(&v)));
    while cc.len() > 4 {
        cc.remove(0);
    }
    Ok(v)
}

fn query_param(q: &str, key: &str) -> Option<String> {
    q.split('&')
        .filter_map(|kv| kv.split_once('='))
        .find(|(k, _)| *k == key)
        .map(|(_, v)| v.to_string())
}

/// A picture of a finished run from any angle:
/// `/api/view?job=..&kind=overlay|distance|plddt&yaw=..&pitch=..&roll=..&zoom=..`. The viewer gets a square frame of
/// fixed size; with `dl` the picture is cropped like the other figures.
fn view_response(st: &State, q: &str) -> Response<Cursor<Vec<u8>>> {
    let id = query_param(q, "job").unwrap_or_default();
    let cp = match comparison_of(st, &id) {
        Ok(v) => v,
        Err(e) => return err(404, &e),
    };
    let f = |k: &str, d: f32| {
        query_param(q, k)
            .and_then(|v| v.parse::<f32>().ok())
            .unwrap_or(d)
    };
    let mut style = cp.1.style.clone();
    style.yaw = f("yaw", style.yaw);
    style.pitch = f("pitch", style.pitch);
    style.roll = f("roll", style.roll);
    style.zoom = f("zoom", style.zoom);
    style.steady = query_param(q, "dl").is_none();
    let kind = match query_param(q, "kind").as_deref() {
        Some("distance") => Kind::Distance,
        Some("plddt") => Kind::Plddt,
        _ => Kind::Overlay,
    };
    let dpi = if query_param(q, "dl").is_some() {
        PRINT_DPI
    } else if query_param(q, "fast").is_some() {
        PREVIEW_DPI / 3.0
    } else {
        PREVIEW_DPI
    };
    let img = figure(&cp.0, kind, &style, dpi);
    let mut r = Response::from_data(encode_jpeg(&img, dpi))
        .with_status_code(200)
        .with_header(Header::from_bytes("Content-Type", "image/jpeg").unwrap());
    if query_param(q, "dl").is_some() {
        let name = match kind {
            Kind::Overlay => "overlay_this_view.jpeg",
            Kind::Distance => "distance_map_this_view.jpeg",
            Kind::Plddt => "confidence_map_this_view.jpeg",
        };
        r = r.with_header(
            Header::from_bytes(
                "Content-Disposition",
                format!("attachment; filename=\"{name}\""),
            )
            .unwrap(),
        );
    }
    r
}

/// Coordinates, colors and the camera of a run, for the page's own fast rotation.
fn viewdata_response(st: &State, q: &str) -> Response<Cursor<Vec<u8>>> {
    let id = query_param(q, "job").unwrap_or_default();
    let cp = match comparison_of(st, &id) {
        Ok(v) => v,
        Err(e) => return err(404, &e),
    };
    let kind = match query_param(q, "kind").as_deref() {
        Some("distance") => Kind::Distance,
        Some("plddt") => Kind::Plddt,
        _ => Kind::Overlay,
    };
    let d = view_data(&cp.0, kind, &cp.1.style);
    let r1 = |x: f32| (x * 100.0).round() / 100.0;
    let chains: Vec<Value> = d
        .chains
        .iter()
        .map(|c| {
            json!({
                "ca": c.ca.iter().flat_map(|p| [r1(p[0]), r1(p[1]), r1(p[2])]).collect::<Vec<_>>(),
                "col": c.col.iter().map(|k| format!("#{:02x}{:02x}{:02x}", k[0], k[1], k[2])).collect::<Vec<_>>(),
                "alpha": c.alpha, "breaks": c.breaks,
            })
        })
        .collect();
    let body = json!({"axes": d.axes, "center": d.center, "half": d.half, "chains": chains});
    Response::from_data(body.to_string())
        .with_status_code(200)
        .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
}

/// Earlier runs, newest first, for the History list.
fn history(st: &State) -> Value {
    let mut rows: Vec<(SystemTime, Value)> = Vec::new();
    if let Ok(rd) = std::fs::read_dir(st.dir()) {
        for e in rd.flatten() {
            let id = e.file_name().to_string_lossy().to_string();
            if !valid_id(&id) {
                continue;
            }
            let p = e.path();
            let Some(sj) = std::fs::read_to_string(p.join("summary.json"))
                .ok()
                .and_then(|t| serde_json::from_str::<Value>(&t).ok())
            else {
                continue;
            };
            let t = std::fs::metadata(p.join("summary.json"))
                .and_then(|m| m.modified())
                .unwrap_or(UNIX_EPOCH);
            rows.push((
                t,
                json!({
                    "id": id, "reference": sj["reference"], "case": sj["case"],
                    "rmsd_all": sj["stats"]["rmsd_all_aligned"], "lddt": sj["stats"]["lddt_ca"],
                    "time": t.duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0),
                }),
            ));
        }
    }
    rows.sort_by_key(|r| std::cmp::Reverse(r.0));
    Value::Array(rows.into_iter().take(50).map(|r| r.1).collect())
}

fn outputs(dir: &Path, c: &Comparison, p: &Params) -> Result<(), String> {
    // one stand-alone 88 mm JPEG per figure at PRINT_DPI; the figures are drawn at the same time
    // (a frame needs about 100 MB while it is drawn), each on its own thread
    let mut kinds = vec![("overlay", Kind::Overlay), ("distance_map", Kind::Distance)];
    if p.style.plddt_panel && c.ref_model.has_conf && c.case_model.has_conf {
        kinds.push(("confidence_map", Kind::Plddt));
    }
    let results: Vec<Result<(), String>> = std::thread::scope(|sc| {
        let mut hs: Vec<_> = kinds
            .iter()
            .map(|&(n, k)| {
                sc.spawn(move || {
                    let img = figure(c, k, &p.style, PRINT_DPI);
                    write(dir, &format!("{n}.jpeg"), &encode_jpeg(&img, PRINT_DPI))
                })
            })
            .collect();
        hs.push(sc.spawn(move || {
            write(
                dir,
                "distance_profile.jpeg",
                &encode_jpeg(&profile_figure(c, PRINT_DPI), PRINT_DPI),
            )
        }));
        hs.into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Err("Internal error while drawing the figure.".into()))
            })
            .collect()
    });
    results.into_iter().collect::<Result<(), String>>()?;
    write(
        dir,
        "residue_distances.csv",
        distance_table_csv(c).as_bytes(),
    )?;
    let mut rob =
        String::from("comparison,rmsd_all_aligned,fit_rmsd,median_distance,max_distance\n");
    for r in &c.robustness {
        rob.push_str(&format!(
            "{},{:.2},{:.2},{:.2},{:.2}\n",
            r.comparison, r.rmsd_all_aligned, r.fit_rmsd, r.median_distance, r.max_distance
        ));
    }
    write(dir, "model_robustness.csv", rob.as_bytes())?;
    write(
        dir,
        if c.case_model.pdb {
            "case_aligned_to_reference.pdb"
        } else {
            "case_aligned_to_reference.cif"
        },
        if c.case_model.pdb {
            transform_pdb(&c.case_model.cif, &c.r, &c.t)
        } else {
            transform_cif(&c.case_model.cif, &c.r, &c.t)
        }
        .as_bytes(),
    )
}

fn valid_id(s: &str) -> bool {
    s.len() == 10 && s.bytes().all(|b| b.is_ascii_hexdigit())
}

fn save_side(
    dir: &Path,
    side: &str,
    data: &[u8],
    files: &HashMap<String, Vec<u8>>,
    fields: &HashMap<String, String>,
) -> Result<(), String> {
    write(dir, &format!("{side}.zip"), data)?;
    write(
        dir,
        &format!("{side}.name"),
        fields
            .get(&format!("{side}_name"))
            .map_or("", |s| s.as_str())
            .as_bytes(),
    )?;
    match files.get(&format!("{side}_pae")).filter(|b| !b.is_empty()) {
        Some(b) => write(dir, &format!("{side}.pae"), b),
        None => Ok(()),
    }
}

fn compare_start(st: &Arc<State>, body: &[u8], ctype: &str) -> Response<Cursor<Vec<u8>>> {
    let parts = multipart(body, ctype);
    let mut fields: HashMap<String, String> = HashMap::new();
    let mut files: HashMap<String, Vec<u8>> = HashMap::new();
    for p in parts {
        if ["ref_zip", "case_zip", "ref_pae", "case_pae"].contains(&p.name.as_str()) {
            files.insert(p.name, p.data);
        } else {
            fields.insert(p.name, String::from_utf8_lossy(&p.data).to_string());
        }
    }
    if let Some(cfg) = &st.server {
        if st.busy.load(Ordering::SeqCst) >= cfg.max_jobs {
            return err(
                503,
                "The server is busy with other comparisons. Try again in a minute.",
            );
        }
        if body.len() > cfg.max_upload {
            return err(413, "The upload is too large for this server.");
        }
    }
    let params = match parse_params(&fields) {
        Ok(p) => p,
        Err(e) => return err(400, &e),
    };
    let n = st.counter.fetch_add(1, Ordering::SeqCst);
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let id = if st.server.is_some() {
        random_hex(5)
    } else {
        format!(
            "{:010x}",
            (stamp as u64 ^ (n as u64).wrapping_mul(0x9E37_79B9)) & 0xFF_FFFF_FFFF
        )
    };
    let dir = st.dir().join(&id);
    if let Err(e) = std::fs::create_dir_all(&dir) {
        return err(500, &format!("Could not create the output folder: {e}"));
    }
    match fields.get("from_job").filter(|s| !s.is_empty()) {
        Some(prev) => {
            if !valid_id(prev) || !st.dir().join(prev).join("ref.zip").exists() {
                return err(
                    400,
                    "The earlier run is no longer available. Upload the files again.",
                );
            }
            for f in SIDE_FILES {
                let src = st.dir().join(prev).join(f);
                if src.exists() {
                    if let Err(e) = std::fs::copy(src, dir.join(f)) {
                        return err(500, &e.to_string());
                    }
                }
            }
        }
        None if fields.get("ref_from").is_some_and(|s| !s.is_empty()) => {
            let prev = fields.get("ref_from").unwrap();
            if !valid_id(prev) || !st.dir().join(prev).join("ref.zip").exists() {
                return err(
                    400,
                    "The earlier run is no longer available. Upload the files again.",
                );
            }
            let Some(c) = files.get("case_zip").filter(|c| !c.is_empty()) else {
                return err(400, "Choose the case file.");
            };
            for f in ["ref.zip", "ref.name", "ref.pae"] {
                let src = st.dir().join(prev).join(f);
                if src.exists() && std::fs::copy(src, dir.join(f)).is_err() {
                    return err(500, "Could not save the uploaded files.");
                }
            }
            if let Err(e) = save_side(&dir, "case", c, &files, &fields) {
                return err(500, &e);
            }
        }
        None => {
            let (Some(r), Some(c)) = (files.get("ref_zip"), files.get("case_zip")) else {
                return err(400, "Choose both files.");
            };
            if r.is_empty() || c.is_empty() {
                return err(400, "Choose both files.");
            }
            if let Err(e) = save_side(&dir, "ref", r, &files, &fields)
                .and(save_side(&dir, "case", c, &files, &fields))
            {
                return err(500, &e);
            }
        }
    }
    st.jobs
        .lock()
        .unwrap()
        .insert(id.clone(), json!({ "status": "running" }));
    st.busy.fetch_add(1, Ordering::SeqCst);
    let (st2, id2) = (Arc::clone(st), id.clone());
    std::thread::spawn(move || {
        let id3 = id2.clone();
        let st3 = Arc::clone(&st2);
        let res = std::panic::catch_unwind(move || run_job(&st3, &id3, &dir, params))
            .unwrap_or_else(|_| Err("Internal error while drawing the figure.".into()));
        let v = match res {
            Ok(s) => json!({ "status": "done", "summary": s }),
            Err(e) => json!({ "status": "error", "error": e }),
        };
        st2.jobs.lock().unwrap().insert(id2, v);
        st2.busy.fetch_sub(1, Ordering::SeqCst);
    });
    json_resp(200, json!({ "job": id }))
}

fn mime(name: &str) -> &'static str {
    match name.rsplit('.').next() {
        Some("png") => "image/png",
        Some("jpeg") | Some("jpg") => "image/jpeg",
        Some("csv") => "text/csv; charset=utf-8",
        Some("json") => "application/json",
        Some("zip") => "application/zip",
        Some("txt") => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}

/// Forget everything of this session: the saved runs (history), the finished jobs and the cached comparisons.
/// The transcript cache stays (it only saves a download). Refused while a comparison is running.
fn reset(st: &State) -> Response<Cursor<Vec<u8>>> {
    if st.busy.load(Ordering::SeqCst) > 0 {
        return err(
            409,
            "A comparison is still running. Wait for it to finish, then reset.",
        );
    }
    let dir = st.dir();
    let ids: Vec<String> = std::fs::read_dir(&dir)
        .map(|rd| {
            rd.flatten()
                .map(|e| e.file_name().to_string_lossy().to_string())
                .filter(|n| valid_id(n))
                .collect()
        })
        .unwrap_or_default();
    if st.server.is_some() {
        let _ = std::fs::remove_dir_all(&dir);
    } else if let Ok(rd) = std::fs::read_dir(&dir) {
        for e in rd.flatten() {
            let _ = if e.path().is_dir() {
                std::fs::remove_dir_all(e.path())
            } else {
                std::fs::remove_file(e.path())
            };
        }
    }
    st.jobs.lock().unwrap().retain(|k, _| !ids.contains(k));
    st.cache_cmp
        .lock()
        .unwrap()
        .retain(|(k, _)| !ids.contains(k));
    json_resp(200, json!({ "cleared": ids.len() }))
}

fn hdr(req: &Request, n: &str) -> String {
    req.headers()
        .iter()
        .find(|h| h.field.as_str().as_str().eq_ignore_ascii_case(n))
        .map(|h| h.value.to_string())
        .unwrap_or_default()
}

fn handle(st: &Arc<State>, mut req: Request) {
    let full = req.url().to_string();
    let (url, query) = full.split_once('?').unwrap_or((full.as_str(), ""));
    let (url, query) = (url.to_string(), query.to_string());
    let method = req.method().clone();
    let ctype = hdr(&req, "Content-Type");
    let mut new_cookie = None;
    if let Some(cfg) = &st.server {
        if url == "/healthz" {
            let _ = req.respond(Response::from_string("ok"));
            return;
        }
        if let Some(pw) = &cfg.password {
            let auth = hdr(&req, "Authorization");
            let given = auth
                .strip_prefix("Basic ")
                .map(|b| b64_decode(b.trim()))
                .unwrap_or_default();
            let given_pw = given
                .split(|&c| c == b':')
                .skip(1)
                .collect::<Vec<_>>()
                .join(&b':');
            if !same_secret(&given_pw, pw.as_bytes()) {
                std::thread::sleep(Duration::from_secs(1));
                let _ = req.respond(
                    Response::from_string("Password needed.")
                        .with_status_code(401)
                        .with_header(
                            Header::from_bytes("WWW-Authenticate", "Basic realm=\"FoldShift\"")
                                .unwrap(),
                        ),
                );
                return;
            }
        }
        let sid = hdr(&req, "Cookie")
            .split(';')
            .filter_map(|c| c.trim().strip_prefix("psid="))
            .find(|v| v.len() == 32 && v.bytes().all(|b| b.is_ascii_hexdigit()))
            .map(str::to_string);
        let sid = sid.unwrap_or_else(|| {
            let v = random_hex(16);
            new_cookie = Some(v.clone());
            v
        });
        SID.with(|s| *s.borrow_mut() = sid);
        if req.body_length().is_some_and(|n| n > cfg.max_upload) {
            let _ = req.respond(err(413, "The upload is too large for this server."));
            return;
        }
    }
    let mut body = Vec::new();
    if method == Method::Post {
        let max = st.server.as_ref().map_or(u64::MAX, |c| c.max_upload as u64);
        let _ = req
            .as_reader()
            .take(max.saturating_add(1))
            .read_to_end(&mut body);
    }
    let now = now_ms(st);
    let resp: Response<Cursor<Vec<u8>>> = match (method, url.as_str()) {
        (Method::Get, "/") => Response::from_string(INDEX)
            .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap()),
        (Method::Get, "/api/version") => json_resp(200, json!({ "version": VERSION })),
        (Method::Post, "/api/ping") => {
            st.last_ping.store(now, Ordering::SeqCst);
            st.bye_at.store(0, Ordering::SeqCst);
            st.seen.store(true, Ordering::SeqCst);
            Response::from_data(Vec::new()).with_status_code(204)
        }
        (Method::Post, "/api/bye") => {
            st.bye_at.store(now + 5000, Ordering::SeqCst);
            Response::from_data(Vec::new()).with_status_code(204)
        }
        (Method::Post, "/api/reset") => reset(st),
        (Method::Get, "/api/history") => json_resp(200, history(st)),
        (Method::Get, "/api/viewdata") => viewdata_response(st, &query),
        (Method::Get, "/api/view") => view_response(st, &query),
        (Method::Post, "/api/sequence") => sequence(st, &body),
        (Method::Post, "/api/compare") => compare_start(st, &body, &ctype),
        (Method::Get, u) if u.starts_with("/api/job/") => {
            let id = &u["/api/job/".len()..];
            let mine = valid_id(id) && st.dir().join(id).is_dir();
            let mem = if mine {
                st.jobs.lock().unwrap().get(id).cloned()
            } else {
                None
            };
            let disk = || {
                let sj = std::fs::read_to_string(st.dir().join(id).join("summary.json")).ok()?;
                let p = std::fs::read_to_string(st.dir().join(id).join("params.json"))
                    .unwrap_or_default();
                Some(
                    json!({ "status": "done", "summary": serde_json::from_str::<Value>(&sj).ok()?,
                    "params": serde_json::from_str::<Value>(&p).unwrap_or(Value::Null) }),
                )
            };
            // a finished run is read from its folder, which also holds the options the page restores
            let unfinished = mem.filter(|v| v["status"] != "done");
            match unfinished.or_else(|| if mine { disk() } else { None }) {
                Some(v) => json_resp(200, v),
                None => err(404, "Unknown job."),
            }
        }
        (Method::Get, u) if u.starts_with("/runs/") => {
            let rest = &u["/runs/".len()..];
            let (id, name) = rest.split_once('/').unwrap_or(("", ""));
            if valid_id(id) && !name.contains('/') && !name.contains("..") && !name.contains('\\') {
                match std::fs::read(st.dir().join(id).join(name)) {
                    Ok(d) => Response::from_data(d)
                        .with_header(Header::from_bytes("Content-Type", mime(name)).unwrap()),
                    Err(_) => Response::from_data(Vec::new()).with_status_code(404),
                }
            } else {
                Response::from_data(Vec::new()).with_status_code(404)
            }
        }
        _ => Response::from_data(Vec::new()).with_status_code(404),
    };
    let mut resp = resp.with_header(Header::from_bytes("Cache-Control", "no-store").unwrap());
    if let (Some(c), Some(cfg)) = (new_cookie, &st.server) {
        resp = resp.with_header(
            Header::from_bytes(
                "Set-Cookie",
                format!(
                    "psid={c}; Path=/; HttpOnly; SameSite=Lax; Max-Age={}",
                    cfg.keep_days * 86400
                ),
            )
            .unwrap(),
        );
    }
    let _ = req.respond(resp);
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--version") {
        println!("FoldShift {VERSION}");
        return;
    }
    let no_browser = args.iter().any(|a| a == "--no-browser");
    let env_num = |k: &str, d: u64| {
        std::env::var(k)
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(d)
    };
    let server_cfg = (args.iter().any(|a| a == "--server")
        || std::env::var("PSPLOT_SERVER").is_ok_and(|v| v == "1"))
    .then(|| ServerCfg {
        password: std::env::var("PSPLOT_PASSWORD")
            .ok()
            .filter(|v| !v.is_empty()),
        max_upload: (env_num("PSPLOT_MAX_UPLOAD_MB", 200) as usize) << 20,
        max_jobs: env_num("PSPLOT_MAX_JOBS", 4).max(1) as usize,
        keep_days: env_num("PSPLOT_KEEP_DAYS", 7).max(1),
    });
    let port: u16 = args
        .iter()
        .position(|a| a == "--port")
        .and_then(|i| args.get(i + 1))
        .and_then(|p| p.parse().ok())
        .unwrap_or(0);
    let home = data_home();
    // Comparisons live in a folder of this session only and are deleted when the app quits.
    // Folders left by a crashed session (older than 12 h) are removed at start; a session that
    // is still open in another window is never touched.
    let sessions = home.join("runs");
    if server_cfg.is_none() {
        remove_old_sessions(&sessions);
    }
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_millis());
    let st = Arc::new(State {
        runs: if server_cfg.is_some() {
            sessions.join("server")
        } else {
            sessions.join(format!("s{}_{stamp}", std::process::id()))
        },
        server: server_cfg,
        cache: home.join("cache"),
        jobs: Mutex::new(HashMap::new()),
        busy: AtomicUsize::new(0),
        last_ping: AtomicI64::new(0),
        bye_at: AtomicI64::new(0),
        seen: AtomicBool::new(false),
        start: Instant::now(),
        counter: AtomicUsize::new(0),
        cache_cmp: Mutex::new(Vec::new()),
    });
    let _ = std::fs::create_dir_all(&st.runs);
    let bind = if st.server.is_some() {
        std::env::var("PSPLOT_BIND")
            .unwrap_or_else(|_| format!("0.0.0.0:{}", if port > 0 { port } else { 8080 }))
    } else {
        format!("127.0.0.1:{port}")
    };
    let server = Arc::new(Server::http(bind).unwrap_or_else(|e| {
        eprintln!("Cannot start the local server: {e}");
        std::process::exit(1);
    }));
    let addr = server
        .server_addr()
        .to_ip()
        .map(|a| a.to_string())
        .unwrap_or_default();
    let url = format!("http://{addr}/");
    println!("FoldShift {VERSION} is running at {url}\nClose the browser tab to quit.");
    if let Some(cfg) = &st.server {
        println!(
            "Server mode: runs are kept {} days per browser session, uploads up to {} MB, {} jobs at once, {}.",
            cfg.keep_days,
            cfg.max_upload >> 20,
            cfg.max_jobs,
            if cfg.password.is_some() { "password required" } else { "no password" }
        );
        let st2 = Arc::clone(&st);
        std::thread::spawn(move || loop {
            expire_sessions(&st2);
            std::thread::sleep(Duration::from_secs(3600));
        });
    }
    for _ in 0..4 {
        let (srv, st) = (Arc::clone(&server), Arc::clone(&st));
        std::thread::spawn(move || {
            while let Ok(req) = srv.recv() {
                let st = Arc::clone(&st);
                std::thread::spawn(move || handle(&st, req));
            }
        });
    }
    if st.server.is_some() {
        loop {
            std::thread::sleep(Duration::from_secs(3600));
        }
    }
    let want_window = cfg!(all(feature = "window", any(windows, target_os = "macos")))
        && !args.iter().any(|a| a == "--browser")
        && !no_browser;
    if want_window {
        let st2 = Arc::clone(&st);
        std::thread::spawn(move || {
            watchdog(&st2);
            remove_session(&st2);
            std::process::exit(0);
        });
        #[cfg(all(feature = "window", any(windows, target_os = "macos")))]
        if native_window(&url) {
            remove_session(&st);
            return; // window closed
        }
        // the window could not be made (for example no WebView2): use the browser instead
        if open::that(&url).is_err() {
            println!("Open this address in a web browser: {url}");
        }
        loop {
            std::thread::sleep(Duration::from_secs(3600));
        }
    }
    if !no_browser && open::that(&url).is_err() {
        println!("Open this address in a web browser: {url}");
    }
    watchdog(&st);
    remove_session(&st);
}

/// Server mode: delete browser-session folders that were not used for `keep_days` days.
fn expire_sessions(st: &State) {
    let Some(cfg) = &st.server else { return };
    let Ok(rd) = std::fs::read_dir(&st.runs) else {
        return;
    };
    for e in rd.flatten() {
        let old = e
            .metadata()
            .and_then(|m| m.modified())
            .ok()
            .and_then(|t| t.elapsed().ok())
            .is_some_and(|d| d > Duration::from_secs(cfg.keep_days * 86400));
        let name = e.file_name().to_string_lossy().to_string();
        if old && name.len() == 32 && e.path().is_dir() && st.busy.load(Ordering::SeqCst) == 0 {
            let _ = std::fs::remove_dir_all(e.path());
        }
    }
}

/// Delete this session's comparisons (history) when the app quits.
fn remove_session(st: &State) {
    if st
        .runs
        .file_name()
        .is_some_and(|n| n.to_string_lossy().starts_with('s'))
    {
        let _ = std::fs::remove_dir_all(&st.runs);
    }
}

/// Remove session folders of earlier crashed sessions (older than 12 h).
fn remove_old_sessions(dir: &Path) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    for e in rd.flatten() {
        let old = e
            .metadata()
            .and_then(|m| m.modified())
            .ok()
            .and_then(|t| t.elapsed().ok())
            .is_some_and(|d| d > Duration::from_secs(12 * 3600));
        if old && e.path().is_dir() {
            let _ = std::fs::remove_dir_all(e.path());
        }
    }
}

/// Quit when the page is closed (never while a figure is being made).
fn watchdog(st: &State) {
    loop {
        std::thread::sleep(Duration::from_millis(500));
        if st.busy.load(Ordering::SeqCst) > 0 {
            continue;
        }
        let now = now_ms(st);
        let bye = st.bye_at.load(Ordering::SeqCst);
        let seen = st.seen.load(Ordering::SeqCst);
        if (bye > 0 && now > bye)
            || (!seen && now > 60_000)
            || (seen && now - st.last_ping.load(Ordering::SeqCst) > 120_000)
        {
            break;
        }
    }
}

/// Show the page in its own window. Returns false when the window could not be created; true when it was closed.
#[cfg(all(feature = "window", any(windows, target_os = "macos")))]
fn native_window(url: &str) -> bool {
    use tao::dpi::LogicalSize;
    use tao::event::{Event, WindowEvent};
    use tao::event_loop::{ControlFlow, EventLoopBuilder};
    use tao::window::WindowBuilder;
    let event_loop = EventLoopBuilder::new().build();
    let Ok(window) = WindowBuilder::new()
        .with_title(format!("FoldShift {VERSION}"))
        .with_inner_size(LogicalSize::new(1100.0, 850.0))
        .build(&event_loop)
    else {
        return false;
    };
    let Ok(_webview) = wry::WebViewBuilder::new().with_url(url).build(&window) else {
        return false;
    };
    event_loop.run(move |event, _, flow| {
        *flow = ControlFlow::Wait;
        if let Event::WindowEvent {
            event: WindowEvent::CloseRequested,
            ..
        } = event
        {
            *flow = ControlFlow::Exit;
        }
        let _ = &_webview;
    });
}
