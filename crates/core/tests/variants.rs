use psplot_core::transcript::{from_sequence, parse_genbank};
use psplot_core::variants::{apply_variants, Transcript};

fn tx() -> Transcript {
    parse_genbank(include_str!("../../../test/data/NM_000558.5.gb")).unwrap()
}
fn v(s: &[&str]) -> Vec<String> {
    s.iter().map(|x| x.to_string()).collect()
}

#[test]
fn reference_protein() {
    let t = tx();
    assert_eq!(t.accession, "NM_000558.5");
    assert_eq!(t.gene, "HBA1");
    // the UniProt P69905 sequence, which RefSeq and Ensembl (ENST00000320868.9) translate to as well
    assert_eq!(t.protein, "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR");
    assert_eq!(
        psplot_core::variants::translate(&t.mrna, t.cds_start).0,
        t.protein
    );
}

#[test]
fn c_and_p_deletion_agree() {
    let t = tx();
    let a = apply_variants(&t, &v(&["c.166_186del"])).unwrap().protein;
    let b = apply_variants(&t, &v(&["p.Val56_Lys62del"]))
        .unwrap()
        .protein;
    assert_eq!(a, b);
    assert_eq!(a.len(), 135);
}

#[test]
fn start_loss_uses_next_met() {
    let r = apply_variants(&tx(), &v(&["c.1A>G"])).unwrap();
    assert!(r.start_loss && r.protein.starts_with("MFLSF") && r.protein.len() == 110);
}

#[test]
fn missense_and_errors() {
    let t = tx();
    // the example of the page: c.179G>A changes codon 60, GGC (Gly) to GAC (Asp)
    let r = apply_variants(&t, &v(&["c.179G>A"])).unwrap();
    assert_eq!(r.protein.as_bytes()[59], b'D');
    let r = apply_variants(&t, &v(&["p.Arg32Gln"])).unwrap();
    assert_eq!(r.protein.as_bytes()[31], b'Q');
    assert!(apply_variants(&t, &v(&["c.5A>G"]))
        .unwrap_err()
        .0
        .contains("reference base"));
    assert!(apply_variants(&t, &v(&["p.Leu113fs"]))
        .unwrap_err()
        .0
        .contains("c. variant"));
    assert!(apply_variants(&t, &v(&["c.1A>G", "p.Arg32Gln"]))
        .unwrap_err()
        .0
        .contains("mix"));
}

#[test]
fn user_sequences() {
    let t = tx();
    let dna: Vec<u8> = t.mrna[t.cds_start..t.cds_end].to_vec();
    let own = from_sequence(&format!(">x\n{}", String::from_utf8(dna).unwrap())).unwrap();
    assert_eq!(own.protein, t.protein);
    assert_eq!(
        apply_variants(&own, &v(&["c.166_186del"])).unwrap().protein,
        apply_variants(&t, &v(&["c.166_186del"])).unwrap().protein
    );
    let prot = from_sequence(&t.protein).unwrap();
    assert_eq!(
        apply_variants(&prot, &v(&["p.Arg32Gln"]))
            .unwrap()
            .protein
            .as_bytes()[31],
        b'Q'
    );
    assert!(apply_variants(&prot, &v(&["c.1A>G"]))
        .unwrap_err()
        .0
        .contains("protein sequence"));
}
