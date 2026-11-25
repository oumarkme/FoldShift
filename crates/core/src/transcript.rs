//! RefSeq transcripts from NCBI (GenBank flat file, cached on disk) and user-pasted sequences.

use crate::variants::{translate, Transcript, VariantError};
use std::path::{Path, PathBuf};

fn err<T>(m: impl Into<String>) -> Result<T, VariantError> {
    Err(VariantError(m.into()))
}

/// Parse a GenBank flat file with one CDS (location `a..b`, no join).
pub fn parse_genbank(text: &str) -> Result<Transcript, VariantError> {
    let mut version = String::new();
    let mut definition = String::new();
    let mut gene = String::new();
    let mut cds: Option<(usize, usize)> = None;
    let mut translation = String::new();
    let mut seq = Vec::new();
    let mut section = "";
    let mut in_def = false;
    let mut in_cds = false;
    let mut in_trans = false;
    for line in text.lines() {
        if let Some(r) = line.strip_prefix("VERSION") {
            version = r.split_whitespace().next().unwrap_or("").to_string();
        }
        if let Some(r) = line.strip_prefix("DEFINITION") {
            definition = r.trim().to_string();
            in_def = true;
            continue;
        }
        if in_def {
            if line.starts_with(' ') {
                definition.push(' ');
                definition.push_str(line.trim());
                continue;
            }
            in_def = false;
        }
        if line.starts_with("FEATURES") {
            section = "features";
            continue;
        }
        if line.starts_with("ORIGIN") {
            section = "origin";
            continue;
        }
        if line.starts_with("//") {
            break;
        }
        match section {
            "features" => {
                if line.len() > 5 && !line.starts_with("      ") && line.starts_with("     ") {
                    // a new feature key
                    in_cds = line[5..].starts_with("CDS");
                    in_trans = false;
                    if in_cds && cds.is_none() {
                        let loc = line[5..].trim_start_matches("CDS").trim();
                        if loc.contains("join") || loc.contains("complement") {
                            return err("This transcript has a complex CDS location (join or complement), which is not supported.");
                        }
                        let loc = loc.replace(['<', '>'], "");
                        let mut it = loc.split("..");
                        let a: usize = it.next().and_then(|x| x.parse().ok()).unwrap_or(0);
                        let b: usize = it.next().and_then(|x| x.parse().ok()).unwrap_or(0);
                        if a == 0 || b < a {
                            return err("Could not read the CDS location.");
                        }
                        cds = Some((a - 1, b));
                    } else {
                        in_cds = false;
                    }
                } else if in_cds {
                    let t = line.trim();
                    if let Some(v) = t.strip_prefix("/gene=") {
                        if gene.is_empty() {
                            gene = v.trim_matches('"').to_string();
                        }
                    } else if let Some(v) = t.strip_prefix("/translation=") {
                        translation = v.trim_matches('"').to_string();
                        in_trans = !v.ends_with('"') || v.len() == 1;
                    } else if in_trans {
                        translation.push_str(t.trim_end_matches('"'));
                        if t.ends_with('"') {
                            in_trans = false;
                        }
                    }
                }
            }
            "origin" => {
                seq.extend(
                    line.bytes()
                        .filter(|b| b.is_ascii_alphabetic())
                        .map(|b| b.to_ascii_uppercase()),
                );
            }
            _ => {}
        }
    }
    let Some((s, e)) = cds else {
        return err("This transcript has no CDS (non-coding RNA).");
    };
    if e > seq.len() {
        return err("The CDS lies outside the sequence in the downloaded record.");
    }
    Ok(Transcript {
        accession: version,
        gene,
        definition,
        mrna: seq,
        cds_start: s,
        cds_end: e,
        protein: translation,
    })
}

pub fn is_accession(s: &str) -> bool {
    let (head, ver) = match s.split_once('.') {
        Some((h, v)) => (h, Some(v)),
        None => (s, None),
    };
    let ok_ver = ver.is_none_or(|v| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()));
    let prefix = ["NM_", "NR_", "XM_", "XR_"]
        .iter()
        .find(|p| head.starts_with(**p));
    prefix.is_some_and(|p| {
        let d = &head[p.len()..];
        !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit())
    }) && ok_ver
}

fn agent(connect: u64, total: u64) -> ureq::Agent {
    use std::time::Duration;
    let b = ureq::Agent::config_builder()
        .timeout_connect(Some(Duration::from_secs(connect)))
        .timeout_global(Some(Duration::from_secs(total)));
    // Windows and macOS: the operating system's TLS (only that provider is compiled in; the
    // default provider would panic). Linux (development): rustls.
    #[cfg(any(windows, target_os = "macos"))]
    let b = b.tls_config(
        ureq::tls::TlsConfig::builder()
            .provider(ureq::tls::TlsProvider::NativeTls)
            .build(),
    );
    b.build().into()
}

/// First source: the NCBI Datasets API (CDS range from the transcript report, mRNA from the gene's FASTA).
fn datasets(accession: &str) -> Result<Transcript, String> {
    use std::io::Read;
    let a = agent(6, 15);
    let base = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene";
    let rep: serde_json::Value = serde_json::from_str(
        &a.get(&format!("{base}/accession/{accession}/product_report"))
            .call()
            .map_err(|e| e.to_string())?
            .body_mut()
            .read_to_string()
            .map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    let prod = rep["reports"][0]["product"].clone();
    let gene_id = prod["gene_id"].as_str().ok_or("no gene id")?.to_string();
    let symbol = prod["symbol"].as_str().unwrap_or("").to_string();
    let t = prod["transcripts"]
        .as_array()
        .and_then(|v| {
            v.iter()
                .find(|t| t["accession_version"].as_str() == Some(accession))
        })
        .ok_or("transcript not in the report")?;
    let r = &t["cds"]["range"][0];
    let num = |k: &str| r[k].as_str().and_then(|x| x.parse::<usize>().ok());
    let (b, e) = (num("begin").ok_or("no CDS")?, num("end").ok_or("no CDS")?);
    if b == 0 || e < b {
        return Err("bad CDS range".into());
    }
    let bytes = a
        .get(&format!(
            "{base}/id/{gene_id}/download?include_annotation_type=FASTA_RNA&filename=g.zip"
        ))
        .call()
        .map_err(|e| e.to_string())?
        .body_mut()
        .read_to_vec()
        .map_err(|e| e.to_string())?;
    let mut z = zip::ZipArchive::new(std::io::Cursor::new(bytes)).map_err(|e| e.to_string())?;
    let mut txt = String::new();
    z.by_name("ncbi_dataset/data/rna.fna")
        .map_err(|e| e.to_string())?
        .read_to_string(&mut txt)
        .map_err(|e| e.to_string())?;
    let mut mrna = Vec::new();
    let mut on = false;
    for l in txt.lines() {
        if let Some(h) = l.strip_prefix('>') {
            on = h.split_whitespace().next() == Some(accession);
        } else if on {
            mrna.extend(l.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if mrna.is_empty() || e > mrna.len() {
        return Err("sequence missing or shorter than the CDS".into());
    }
    let (protein, _) = translate(&mrna, b - 1);
    let name = t["name"].as_str().unwrap_or("");
    Ok(Transcript {
        accession: accession.to_string(),
        gene: symbol.clone(),
        definition: format!("{symbol}, {name} ({accession})"),
        mrna,
        cds_start: b - 1,
        cds_end: e,
        protein,
    })
}

fn download(accession: &str) -> Result<String, VariantError> {
    let url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gb&retmode=text&tool=psplot"
    );
    let agent = agent(6, 20);
    let mut last = String::new();
    for _ in 0..1 {
        match agent.get(&url).call() {
            Ok(mut resp) => match resp.body_mut().read_to_string() {
                Ok(t) => return Ok(t),
                Err(e) => last = e.to_string(),
            },
            Err(e) => last = e.to_string(),
        }
    }
    err(last)
}

/// Ensembl transcript ID such as ENST00000320868 or ENST00000320868.7.
pub fn is_ensembl(s: &str) -> bool {
    let h = s.split_once('.').map_or(s, |x| x.0);
    h.len() == 15 && h.starts_with("ENST") && h[4..].bytes().all(|b| b.is_ascii_digit())
}

/// UniProt accession such as P69905 or P69905-2 (isoform).
/// Official pattern: [OPQ][0-9][A-Z0-9]{3}[0-9] or [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}.
pub fn is_uniprot(s: &str) -> bool {
    let (h, iso) = match s.split_once('-') {
        Some((h, i)) => (h, Some(i)),
        None => (s, None),
    };
    if iso.is_some_and(|i| i.is_empty() || !i.bytes().all(|c| c.is_ascii_digit())) {
        return false;
    }
    let b = h.as_bytes();
    let al = |c: u8| c.is_ascii_uppercase() || c.is_ascii_digit();
    let grp = |x: &[u8]| x[0].is_ascii_uppercase() && al(x[1]) && al(x[2]) && x[3].is_ascii_digit();
    if b.len() == 6 && matches!(b[0], b'O' | b'P' | b'Q') {
        return b[1].is_ascii_digit() && b[2..5].iter().all(|&c| al(c)) && b[5].is_ascii_digit();
    }
    (b.len() == 6 || b.len() == 10)
        && b[0].is_ascii_uppercase()
        && !matches!(b[0], b'O' | b'P' | b'Q')
        && b[1].is_ascii_digit()
        && grp(&b[2..6])
        && (b.len() == 6 || grp(&b[6..10]))
}

fn http_text(url: &str) -> Result<String, VariantError> {
    agent(6, 20)
        .get(url)
        .call()
        .map_err(|e| VariantError(format!("Could not reach the server ({e}). Check the internet connection or proxy, or paste the sequence instead.")))?
        .body_mut()
        .read_to_string()
        .map_err(|e| VariantError(format!("Could not read the answer ({e}).")))
}

/// Ensembl REST: cDNA, CDS (located inside the cDNA) and the gene name of one transcript.
fn ensembl(id: &str) -> Result<Transcript, VariantError> {
    let seq = |t: &str| {
        http_text(&format!(
            "https://rest.ensembl.org/sequence/id/{id}?type={t};content-type=text/plain"
        ))
    };
    let mrna: Vec<u8> = seq("cdna")?
        .trim()
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let cds: Vec<u8> = seq("cds")?
        .trim()
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if mrna.len() < 10 || cds.len() < 3 || !mrna.iter().all(u8::is_ascii_uppercase) {
        return err(format!("Ensembl has no coding sequence for {id}."));
    }
    let Some(start) = mrna.windows(cds.len()).position(|w| w == cds.as_slice()) else {
        return err("The coding sequence was not found inside the transcript.");
    };
    let name = serde_json::from_str::<serde_json::Value>(&http_text(&format!(
        "https://rest.ensembl.org/lookup/id/{id}?content-type=application/json"
    ))?)
    .ok()
    .and_then(|v| v["display_name"].as_str().map(String::from))
    .unwrap_or_default();
    let gene = name
        .rsplit_once('-')
        .map_or(name.clone(), |x| x.0.to_string());
    let (protein, _) = translate(&mrna, start);
    Ok(Transcript {
        accession: id.to_string(),
        gene,
        definition: format!("{name} ({id}, Ensembl)"),
        mrna,
        cds_start: start,
        cds_end: start + cds.len(),
        protein,
    })
}

/// UniProt REST: protein only (use p. variants).
fn uniprot(acc: &str) -> Result<Transcript, VariantError> {
    let text = http_text(&format!("https://rest.uniprot.org/uniprotkb/{acc}.fasta"))?;
    let head = text.lines().next().unwrap_or("").to_string();
    let protein: String = text.lines().skip(1).map(str::trim).collect();
    if !head.starts_with('>') || protein.is_empty() {
        return err(format!("UniProt has no entry {acc}."));
    }
    let gene = head
        .split(" GN=")
        .nth(1)
        .and_then(|x| x.split_whitespace().next())
        .unwrap_or("")
        .to_string();
    let title = head.split(" OS=").next().unwrap_or("");
    let title = title.split_once(' ').map_or("", |x| x.1);
    Ok(Transcript {
        accession: acc.to_string(),
        gene,
        definition: format!("{title} (UniProt {acc}; protein only, use p. variants)"),
        mrna: vec![],
        cds_start: 0,
        cds_end: 0,
        protein,
    })
}

pub fn get_transcript(accession: &str, cache: &Path) -> Result<Transcript, VariantError> {
    let acc = accession.trim();
    if is_ensembl(acc) {
        return ensembl(acc);
    }
    if is_uniprot(acc) {
        return uniprot(acc);
    }
    if !is_accession(acc) {
        return err(format!(
            "'{acc}' is not recognized. Use a RefSeq transcript (NM_000558.5), an Ensembl transcript (ENST00000320868.9) or a UniProt accession (P69905)."
        ));
    }
    let gb: PathBuf = cache.join(format!("{acc}.gb"));
    let ds: PathBuf = cache.join(format!("{acc}.ds.tsv"));
    if let Ok(t) = std::fs::read_to_string(&gb) {
        return parse_genbank(&t);
    }
    if let Some(t) = std::fs::read_to_string(&ds)
        .ok()
        .and_then(|t| read_ds_cache(acc, &t))
    {
        return Ok(t);
    }
    let _ = std::fs::create_dir_all(cache);
    // 1. NCBI Datasets API; 2. NCBI E-utilities (GenBank record), started in parallel if the
    // first has not answered after 5 s. The first success wins; the whole lookup ends within 25 s.
    use std::sync::mpsc;
    use std::time::{Duration, Instant};
    enum Got {
        Ds(Result<Transcript, String>),
        Gb(Result<String, VariantError>),
    }
    let (tx, rx) = mpsc::channel();
    let (a1, tx1) = (acc.to_string(), tx.clone());
    std::thread::spawn(move || {
        let _ = tx1.send(Got::Ds(datasets(&a1)));
    });
    let (a2, tx2) = (acc.to_string(), tx);
    let started = Instant::now();
    let mut gb_started = false;
    let mut ds_err = String::new();
    let mut gb_err = String::new();
    let (mut ds_done, mut gb_done) = (false, false);
    let mut tx2 = Some(tx2);
    loop {
        if !gb_started && (ds_done || started.elapsed() > Duration::from_secs(5)) {
            gb_started = true;
            if let Some(t2) = tx2.take() {
                let a = a2.clone();
                std::thread::spawn(move || {
                    let _ = t2.send(Got::Gb(download(&a)));
                });
            }
        }
        if ds_done && gb_done {
            break;
        }
        if started.elapsed() > Duration::from_secs(25) {
            break;
        }
        match rx.recv_timeout(Duration::from_millis(200)) {
            Ok(Got::Ds(Ok(t))) => {
                let line = format!(
                    "{}\t{}\t{}\t{}\n{}\n",
                    t.gene,
                    t.cds_start,
                    t.cds_end,
                    t.definition,
                    String::from_utf8_lossy(&t.mrna)
                );
                let _ = std::fs::write(&ds, line);
                return Ok(t);
            }
            Ok(Got::Ds(Err(e))) => {
                ds_done = true;
                ds_err = e;
            }
            Ok(Got::Gb(Ok(text))) => {
                if let Ok(t) = parse_genbank(&text) {
                    let _ = std::fs::write(&gb, &text);
                    return Ok(t);
                }
                gb_done = true;
                gb_err = "the record could not be read".into();
            }
            Ok(Got::Gb(Err(VariantError(m)))) => {
                gb_done = true;
                gb_err = m;
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                // a lookup thread ended without an answer (for example a panic)
                if !ds_done && ds_err.is_empty() {
                    ds_err = "stopped unexpectedly".into();
                }
                if gb_started && gb_err.is_empty() {
                    gb_err = "stopped unexpectedly".into();
                }
                break;
            }
            Err(_) => {}
        }
    }
    err(format!(
        "Could not get {acc} from NCBI. Datasets API: {}. E-utilities: {}. Check the internet connection or proxy, or paste the sequence instead.",
        if ds_err.is_empty() { "no answer in time" } else { &ds_err },
        if gb_err.is_empty() { "no answer in time" } else { &gb_err }
    ))
}

fn read_ds_cache(acc: &str, text: &str) -> Option<Transcript> {
    let mut l = text.lines();
    let h: Vec<&str> = l.next()?.splitn(4, '\t').collect();
    let mrna: Vec<u8> = l.next()?.bytes().collect();
    let (s, e) = (
        h.get(1)?.parse::<usize>().ok()?,
        h.get(2)?.parse::<usize>().ok()?,
    );
    if e > mrna.len() || s >= e {
        return None;
    }
    Some(Transcript {
        accession: acc.to_string(),
        gene: h[0].to_string(),
        definition: h.get(3)?.to_string(),
        protein: translate(&mrna, s).0,
        mrna,
        cds_start: s,
        cds_end: e,
    })
}

/// Build a Transcript from pasted FASTA or raw text. DNA/RNA: c.1 is the A of the first ATG. Protein: p. variants only.
pub fn from_sequence(text: &str) -> Result<Transcript, VariantError> {
    let mut name = "user_sequence".to_string();
    let mut lines = text.trim().lines().map(str::trim);
    let mut body: Vec<&str> = Vec::new();
    if let Some(first) = lines.next() {
        if let Some(h) = first.strip_prefix('>') {
            if let Some(n) = h.split_whitespace().next() {
                name = n.to_string();
            }
        } else {
            body.push(first);
        }
    }
    body.extend(lines);
    let seq: Vec<u8> = body
        .concat()
        .bytes()
        .filter(|b| !b.is_ascii_whitespace() && !b.is_ascii_digit())
        .map(|b| b.to_ascii_uppercase())
        .map(|b| if b == b'U' { b'T' } else { b })
        .collect();
    if seq.len() < 10 {
        return err("The sequence is too short (at least 10 letters).");
    }
    if !seq.iter().all(|b| b.is_ascii_uppercase() || *b == b'*') {
        return err("The sequence contains characters that are not letters.");
    }
    if seq.iter().all(|b| b"ACGTN".contains(b)) {
        let Some(start) = seq.windows(3).position(|w| w == b"ATG") else {
            return err("No ATG start codon was found in the DNA sequence. Paste it starting at the start codon, or paste the protein sequence.");
        };
        let (protein, stopped) = translate(&seq, start);
        let end = start + (protein.len() + usize::from(stopped)) * 3;
        return Ok(Transcript {
            accession: name.clone(),
            gene: String::new(),
            definition: format!(
                "{name} (DNA you pasted; c.1 is the first ATG, at position {})",
                start + 1
            ),
            mrna: seq,
            cds_start: start,
            cds_end: end,
            protein,
        });
    }
    let protein: String = seq
        .iter()
        .map(|&b| b as char)
        .collect::<String>()
        .trim_end_matches('*')
        .to_string();
    Ok(Transcript {
        accession: name.clone(),
        gene: String::new(),
        definition: format!("{name} (protein you pasted; use p. variants)"),
        mrna: vec![],
        cds_start: 0,
        cds_end: 0,
        protein,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn id_kinds() {
        for u in ["P69905", "P69905-2", "P04637", "A0A024R161", "Q9UKR5"] {
            assert!(is_uniprot(u), "{u}");
        }
        for u in ["NM_000558.5", "ENST00000320868", "HBA1", "O6068", "P69905-"] {
            assert!(!is_uniprot(u), "{u}");
        }
        assert!(is_ensembl("ENST00000320868.7") && !is_ensembl("NM_000558.5"));
    }
}
