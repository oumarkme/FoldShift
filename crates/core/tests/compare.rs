use psplot_core::compare::{compare, Options};
use psplot_core::structure::load_zip;
use std::fs::File;

fn load(n: &str) -> psplot_core::structure::AfResult {
    load_zip(
        File::open(format!(
            "{}/../../test/data/{n}.zip",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn deletion_and_hinge() {
    let (r, f) = (load("hba1_ref"), load("hba1_del_hinge"));
    let c = compare(
        &r,
        &f,
        &Options {
            robustness: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(c.missing_in_case, vec![(56, 62)]);
    let loc: Vec<f64> = c
        .blocks
        .iter()
        .map(|b| (b.rmsd_local_fit * 10.0).round() / 10.0)
        .collect();
    // each rigid half fits itself to about the noise of the test models (0.2 A per coordinate)
    assert_eq!(loc, vec![0.5, 1.0]);
    // the first half sits 25 degrees away from where a global fit on the second half puts it
    assert!(c.blocks[0].rmsd_global_fit > 5.0 && c.blocks[1].rmsd_global_fit < 1.5);
    assert!(c.stats.fit_rmsd < 1.0 && c.stats.fit_residues > 70);
}

#[test]
fn n_terminal_truncation() {
    let (r, m) = (load("hba1_ref"), load("hba1_trunc"));
    let c = compare(
        &r,
        &m,
        &Options {
            robustness: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(c.missing_in_case, vec![(1, 39)]);
    assert_eq!(c.blocks.len(), 1);
    assert!((c.blocks[0].rmsd_local_fit - 0.5).abs() < 0.05);
}

#[test]
fn custom_region_and_robustness() {
    let (r, f) = (load("hba1_ref"), load("hba1_del_hinge"));
    let c = compare(
        &r,
        &f,
        &Options {
            fit_region: "72-142".into(),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(c.stats.fit_rmsd < 1.0);
    assert_eq!(c.robustness.len(), 5 + 4);
    assert!(compare(
        &r,
        &f,
        &Options {
            fit_region: "1-5".into(),
            ..Default::default()
        }
    )
    .is_err());
}
