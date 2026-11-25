//! Predicted aligned error (PAE) matrices from JSON, NumPy .npy or .npz files.
//! Pickle files are refused: opening an untrusted pickle can run code.

use std::io::Read;

/// An n x n matrix in row-major order, in Angstrom.
#[derive(Debug)]
pub struct Pae {
    pub n: usize,
    pub v: Vec<f32>,
}

impl Pae {
    /// The matrix of the residues at `idx` (rows and columns), in that order.
    pub fn subset(&self, idx: &[usize]) -> Pae {
        let n = idx.len();
        let mut v = Vec::with_capacity(n * n);
        for &i in idx {
            for &j in idx {
                v.push(self.v[i * self.n + j]);
            }
        }
        Pae { n, v }
    }

    pub fn mean(&self) -> f64 {
        self.v.iter().map(|&x| x as f64).sum::<f64>() / self.v.len().max(1) as f64
    }
    pub fn max(&self) -> f32 {
        self.v.iter().cloned().fold(0.0, f32::max)
    }
}

fn square(rows: Vec<Vec<f32>>) -> Result<Pae, String> {
    let n = rows.len();
    if n < 2 || rows.iter().any(|r| r.len() != n) {
        return Err("The PAE matrix is not square.".into());
    }
    Ok(Pae {
        n,
        v: rows.into_iter().flatten().collect(),
    })
}

fn matrix(v: &serde_json::Value) -> Option<Vec<Vec<f32>>> {
    v.as_array()?
        .iter()
        .map(|r| {
            r.as_array()?
                .iter()
                .map(|x| x.as_f64().map(|f| f as f32))
                .collect::<Option<Vec<f32>>>()
        })
        .collect()
}

fn from_json(text: &str) -> Result<Pae, String> {
    let v: serde_json::Value =
        serde_json::from_str(text).map_err(|_| "This is not a readable JSON file.".to_string())?;
    // AlphaFold Server and ColabFold: {"pae": [[...]]}; AlphaFold database: [{"predicted_aligned_error": [[...]]}]
    let obj = if v.is_array() && v[0].is_object() {
        &v[0]
    } else {
        &v
    };
    for key in ["pae", "predicted_aligned_error"] {
        if let Some(m) = obj.get(key).and_then(matrix) {
            return square(m);
        }
    }
    if let Some(m) = matrix(obj) {
        return square(m);
    }
    Err("No PAE matrix was found in this JSON file (look for the file with \"pae\" or \"predicted_aligned_error\", for example fold_..._full_data_0.json).".into())
}

fn from_npy(b: &[u8]) -> Result<Pae, String> {
    let bad = || "This .npy file could not be read.".to_string();
    if b.len() < 12 || &b[..6] != b"\x93NUMPY" {
        return Err(bad());
    }
    let (hl, start) = if b[6] == 1 {
        (u16::from_le_bytes([b[8], b[9]]) as usize, 10)
    } else {
        (u32::from_le_bytes([b[8], b[9], b[10], b[11]]) as usize, 12)
    };
    let head = std::str::from_utf8(b.get(start..start + hl).ok_or_else(bad)?).map_err(|_| bad())?;
    let descr = head
        .split("'descr':")
        .nth(1)
        .and_then(|x| x.split('\'').nth(1))
        .ok_or_else(bad)?;
    let fortran = head.contains("'fortran_order': True");
    let shape: Vec<usize> = head
        .split("'shape':")
        .nth(1)
        .and_then(|x| x.split('(').nth(1))
        .and_then(|x| x.split(')').next())
        .ok_or_else(bad)?
        .split(',')
        .filter_map(|t| t.trim().parse().ok())
        .collect();
    let data = &b[start + hl..];
    let vals: Vec<f32> = match descr {
        "<f4" => data
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| f32::from_le_bytes(*c))
            .collect(),
        "<f8" => data
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| f64::from_le_bytes(*c) as f32)
            .collect(),
        _ => {
            return Err(format!(
            "The number type {descr} in this .npy file is not supported (use float32 or float64)."
        ))
        }
    };
    let n = match shape.as_slice() {
        [a, b2] if a == b2 => *a,
        _ => return Err("The .npy array is not a square matrix.".into()),
    };
    if vals.len() != n * n {
        return Err(bad());
    }
    let v = if fortran {
        (0..n * n).map(|k| vals[(k % n) * n + k / n]).collect()
    } else {
        vals
    };
    Ok(Pae { n, v })
}

/// Read a PAE matrix from the bytes of a .json, .npy or .npz file (the kind is detected from the content).
pub fn parse_pae(b: &[u8]) -> Result<Pae, String> {
    if b.starts_with(b"\x93NUMPY") {
        return from_npy(b);
    }
    if b.starts_with(b"PK") {
        let mut z = zip::ZipArchive::new(std::io::Cursor::new(b))
            .map_err(|_| "This archive could not be read.".to_string())?;
        let names: Vec<String> = (0..z.len())
            .filter_map(|i| z.by_index(i).ok().map(|f| f.name().to_string()))
            .collect();
        let pick = names
            .iter()
            .find(|n| n.to_ascii_lowercase().contains("pae") && n.ends_with(".npy"))
            .or_else(|| names.iter().find(|n| n.ends_with(".npy")))
            .ok_or("No PAE array (.npy) was found in this archive.")?;
        let mut buf = Vec::new();
        z.by_name(pick)
            .map_err(|e| e.to_string())?
            .read_to_end(&mut buf)
            .map_err(|e| e.to_string())?;
        return from_npy(&buf);
    }
    if b.first() == Some(&0x80) {
        return Err("Pickle files (.pkl) are not supported because opening one can run code. Save the PAE as .json or .npy instead.".into());
    }
    from_json(std::str::from_utf8(b).map_err(|_| "This file is not a PAE file.".to_string())?)
}
