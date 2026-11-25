//! Tests of the v0.0.2 additions: secondary structure, BLOSUM62 pairing, several chains, substituted residues.

use psplot_core::compare::{align_sequences, compare, substituted_residues, Options};
use psplot_core::ss::{assign_ss, Ss};
use psplot_core::structure::{load_zip, parse_cif};
use std::fs::File;

fn zip(n: &str) -> psplot_core::structure::AfResult {
    load_zip(
        File::open(format!(
            "{}/../../test/data/{n}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn ideal_helix_and_strand() {
    // alpha helix: radius 2.3 A, rise 1.5 A, 100 degrees per residue
    let helix: Vec<[f64; 3]> = (0..14)
        .map(|i| {
            let a = (i as f64 * 100.0).to_radians();
            [2.3 * a.cos(), 2.3 * a.sin(), 1.5 * i as f64]
        })
        .collect();
    let ss = assign_ss(&helix, &vec![false; helix.len()]);
    assert!(ss.iter().filter(|s| **s == Ss::Helix).count() >= 10);
    // extended strand: zig-zag, 3.3 A per residue along x
    let strand: Vec<[f64; 3]> = (0..12)
        .map(|i| [3.3 * i as f64, if i % 2 == 0 { 0.0 } else { 0.9 }, 0.0])
        .collect();
    let ss = assign_ss(&strand, &vec![false; strand.len()]);
    assert!(ss.iter().filter(|s| **s == Ss::Strand).count() >= 8);
}

#[test]
fn blosum_pairs_through_substitutions() {
    let a = b"MKTAYIAKQRQISFVKSHFSRQ";
    let b = b"MKTAYLAKQRQISFVKSHFSRQ"; // I6L: similar residues
    for blosum in [true, false] {
        let p = align_sequences(a, b, blosum);
        assert_eq!(p.len(), a.len());
    }
    assert_eq!(substituted_residues("MKTAYIAKQ", "MKTAYLAKQ"), vec![6]);
    assert!(substituted_residues("MKTAYIAKQ", "MKTAYIAKQ").is_empty());
}

#[test]
fn two_chains_are_parsed_and_never_bridged() {
    let cif = "data_x\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\nATOM 1 CA ALA A 1 0.0 0.0 0.0 90.0\nATOM 2 CA GLY A 2 3.8 0.0 0.0 80.0\nATOM 3 CA SER B 1 4.0 0.0 0.0 70.0\nATOM 4 CA LYS B 2 7.8 0.0 0.0 60.0\n#\n";
    let d = parse_cif(cif).unwrap();
    assert_eq!(d.aa, vec!['A', 'G', 'S', 'K']);
    assert_eq!(d.chain, vec!["A", "A", "B", "B"]);
    // residue 2 (chain A) and 3 (chain B) are 0.2 A apart, but there is no bond between chains
    assert_eq!(d.breaks, vec![false, true, false, false]);
}

#[test]
fn new_statistics_are_sane() {
    let (r, c) = (zip("hba1_ref.zip"), zip("hba1_del_hinge.zip"));
    let cmp = compare(&r, &c, &Options::default()).unwrap();
    let s = &cmp.stats;
    assert!((0.0..=1.0).contains(&s.lddt_ca));
    assert_eq!((s.chains_ref, s.chains_case), (1, 1));
    assert_eq!(cmp.lddt.len(), cmp.pairs.len());
    assert!(cmp.low_conf.iter().filter(|x| **x).count() == s.low_confidence_shifts);
    // the global-fit numbers do not depend on the pairing choice for this case (no substitutions)
    let id = compare(
        &r,
        &c,
        &Options {
            blosum: false,
            ..Options::default()
        },
    )
    .unwrap();
    assert_eq!(id.pairs, cmp.pairs);
}

#[test]
fn pdb_and_pae_files() {
    use psplot_core::structure::load_any;
    let z = zip("hba1_del_hinge.zip");
    let m = &z.models[z.best()];
    // build a PDB file from the C-alpha atoms; B-factor on the 0-1 scale
    let mut pdb = String::new();
    for i in 0..m.ca.len() {
        pdb.push_str(&format!(
            "ATOM  {:5} {:<4} {:3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C\n",
            i + 1,
            " CA ".trim(),
            "ALA",
            m.chain[i],
            m.resnums[i],
            m.ca[i][0],
            m.ca[i][1],
            m.ca[i][2],
            m.plddt[i] / 100.0
        ));
    }
    let r = load_any(pdb.as_bytes(), "case.pdb", None).unwrap();
    assert_eq!(r.models[0].ca.len(), m.ca.len());
    assert!(r.models[0].has_conf);
    assert!((r.models[0].plddt[10] - m.plddt[10]).abs() < 0.6);
    assert_eq!(r.source, "structure file");
    // npy PAE of the right size, and a wrong size refused
    let n = m.ca.len();
    let mut npy = b"\x93NUMPY\x01\x00".to_vec();
    let h = format!("{{'descr': '<f4', 'fortran_order': False, 'shape': ({n}, {n}), }}");
    let pad = 64 - (10 + h.len() + 1) % 64;
    let h = format!("{h}{}\n", " ".repeat(pad));
    npy.extend((h.len() as u16).to_le_bytes());
    npy.extend(h.as_bytes());
    for _ in 0..n * n {
        npy.extend(3.5f32.to_le_bytes());
    }
    let r = load_any(pdb.as_bytes(), "case.pdb", Some(&npy)).unwrap();
    assert!((r.models[0].pae.as_ref().unwrap().mean() - 3.5).abs() < 1e-6);
    assert!(load_any(pdb.as_bytes(), "x.pdb", Some(b"\x80\x04pickle")).is_err());
    let small = br#"{"pae": [[0,1],[1,0]]}"#;
    assert!(load_any(pdb.as_bytes(), "x.pdb", Some(small)).is_err());
    // the AlphaFold Server zip carries its own PAE
    assert!(m.pae.is_some());
}

#[test]
fn experimental_b_factors_are_not_confidence() {
    use psplot_core::structure::{is_experimental, load_any};
    let atoms: String = (1..=12)
        .map(|i| {
            format!(
                "ATOM  {:5}  CA  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C\n",
                i,
                i,
                i as f64 * 3.8,
                0.0,
                0.0,
                20.0 + i as f64 * 3.0
            )
        })
        .collect();
    let xray = format!("HEADER    TEST\nEXPDTA    X-RAY DIFFRACTION\n{atoms}END\n");
    assert!(is_experimental(&xray));
    let r = load_any(xray.as_bytes(), "x.pdb", None).unwrap();
    assert!(!r.models[0].has_conf);
    assert!(r.models[0].plddt.iter().all(|&v| v == 100.0));
    // the same numbers without an experiment line are read as pLDDT
    let pred = format!("HEADER    TEST\n{atoms}END\n");
    assert!(!is_experimental(&pred));
    let r = load_any(pred.as_bytes(), "y.pdb", None).unwrap();
    assert!(r.models[0].has_conf && r.models[0].plddt[0] < 30.0);
    assert!(is_experimental(
        "data_x\n_exptl.method 'ELECTRON MICROSCOPY'\n"
    ));
    assert!(!is_experimental(
        "data_x\n_exptl.method 'THEORETICAL MODEL'\n"
    ));
}

/// A small RCSB-style PDB text: `resnums` with C-alpha on a gently curved line, one chain letter.
fn pdb_text(chain: char, resnums: &[i32], aa: &[&str], shift: f64) -> String {
    let mut s = String::from("HEADER    TEST\nEXPDTA    X-RAY DIFFRACTION\n");
    for (k, (&n, a)) in resnums.iter().zip(aa).enumerate() {
        let t = n as f64 * 1.0;
        s.push_str(&format!(
            "ATOM  {:5}  CA  {} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C\n",
            k + 1,
            a,
            chain,
            n,
            3.8 * t.cos() * 2.0 + shift,
            3.8 * t.sin() * 2.0,
            t * 1.5,
            30.0
        ));
    }
    s + "END\n"
}

#[test]
fn chains_and_numbers() {
    use psplot_core::compare::{compare, Options};
    use psplot_core::structure::load_any;
    // 30 residues numbered 100..129; the case lacks 110-113 (a disordered loop) and has chain B too
    let names = [
        "ALA", "GLY", "LEU", "VAL", "SER", "LYS", "GLU", "ASP", "THR", "PHE",
    ];
    let nums: Vec<i32> = (100..130).collect();
    let aa: Vec<&str> = nums.iter().map(|n| names[(*n as usize * 7) % 10]).collect();
    let r_txt = pdb_text('A', &nums, &aa, 0.0);
    let keep: Vec<usize> = (0..30).filter(|i| !(10..14).contains(i)).collect();
    let c_nums: Vec<i32> = keep.iter().map(|&i| nums[i]).collect();
    let c_aa: Vec<&str> = keep.iter().map(|&i| aa[i]).collect();
    let c_txt = pdb_text('A', &c_nums, &c_aa, 0.2) + &pdb_text('B', &c_nums, &c_aa, 40.0);
    let r = load_any(r_txt.as_bytes(), "r.pdb", None).unwrap();
    let mut c = load_any(c_txt.as_bytes(), "c.pdb", None).unwrap();
    assert_eq!(c.models[0].ca.len(), 52);
    assert!(c.keep_chains("Z").unwrap_err().contains("chains: A, B"));
    c.keep_chains("A").unwrap();
    assert_eq!(c.models[0].ca.len(), 26);
    assert_eq!(c.sequence.len(), 26);
    let cmp = compare(&r, &c, &Options::default()).unwrap();
    assert!(cmp.paired_by_number);
    assert_eq!(cmp.pairs.len(), 26);
    // residue numbers are the numbers of the file, and the gap shows as absent residues 110-113
    assert_eq!(cmp.ref_num[0], 100);
    assert_eq!(cmp.missing_in_case, vec![(11, 14)]);
    assert!(psplot_core::compare::plain_summary(&cmp)
        .iter()
        .any(|l| l.contains("110-113")));
}

#[test]
fn rcsb_style_cif_columns() {
    use psplot_core::structure::parse_cif;
    // trailing blanks after the column names, author numbering, a second alternative location, a second model
    let t = "data_x\nloop_\n_atom_site.group_PDB \n_atom_site.id \n_atom_site.type_symbol \n_atom_site.label_atom_id \n_atom_site.label_alt_id \n_atom_site.label_comp_id \n_atom_site.label_asym_id \n_atom_site.label_entity_id \n_atom_site.label_seq_id \n_atom_site.pdbx_PDB_ins_code \n_atom_site.Cartn_x \n_atom_site.Cartn_y \n_atom_site.Cartn_z \n_atom_site.occupancy \n_atom_site.B_iso_or_equiv \n_atom_site.pdbx_formal_charge \n_atom_site.auth_seq_id \n_atom_site.auth_comp_id \n_atom_site.auth_asym_id \n_atom_site.auth_atom_id \n_atom_site.pdbx_PDB_model_num \nATOM 1 C CA . THR A 1 2 ? 0.0 0.0 0.0 1.0 66.7 ? 389 THR C CA 1\nATOM 2 C CA A SER A 1 3 ? 3.8 0.0 0.0 1.0 60.0 ? 390 SER C CA 1\nATOM 3 C CA B SER A 1 3 ? 3.9 0.0 0.0 1.0 60.0 ? 390 SER C CA 1\nATOM 4 C CA . GLY A 1 4 ? 7.6 0.0 0.0 1.0 60.0 ? 391 GLY C CA 1\nATOM 5 C CA . ALA A 1 5 ? 0.0 0.0 0.0 1.0 60.0 ? 392 ALA C CA 2\n#\n";
    let d = parse_cif(t).unwrap();
    assert_eq!(d.resnums, vec![389, 390, 391]);
    assert_eq!(d.chain, vec!["C", "C", "C"]);
    assert_eq!(d.aa.iter().collect::<String>(), "TSG");
}

#[test]
fn open_pdb_pair_without_confidence() {
    use psplot_core::structure::load_any;
    let read = |n: &str| {
        std::fs::read(format!(
            "{}/../../test/data/{n}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap()
    };
    // T4 lysozyme, two experimental entries (RCSB 1L63 and 1L90): no pLDDT, so every residue counts as confident
    let (r, c) = (
        load_any(&read("1L63.pdb"), "1L63.pdb", None).unwrap(),
        load_any(&read("1L90.pdb"), "1L90.pdb", None).unwrap(),
    );
    assert!(!r.models[0].has_conf && !c.models[0].has_conf);
    let cmp = compare(&r, &c, &Options::default()).unwrap();
    // point mutants of the same fold: almost identical backbones
    assert!(cmp.stats.fit_rmsd < 0.5 && cmp.stats.fit_residues > 150);
}

#[test]
fn alphafold_server_hba1_g60d() {
    // real AlphaFold Server results: HBA1 reference and the case c.179G>A (p.Gly60Asp), five models each
    let (r, c) = (zip("fold_hba_reference.zip"), zip("fold_hba_case.zip"));
    assert_eq!(
        (r.source.as_str(), r.sequence.len(), r.models.len()),
        ("AlphaFold Server", 142, 5)
    );
    assert_eq!(substituted_residues(&r.sequence, &c.sequence), vec![60]);
    assert!(r.models[r.best()].pae.is_some() && c.models[c.best()].pae.is_some());
    let cmp = compare(&r, &c, &Options::default()).unwrap();
    // one substitution in a well-folded helical protein: one block, nearly identical backbones
    assert!(cmp.missing_in_case.is_empty() && cmp.blocks.len() == 1);
    assert!(cmp.stats.fit_rmsd < 0.5 && cmp.stats.fit_residues > 130);
    assert!(cmp.stats.lddt_ca > 0.98);
}
