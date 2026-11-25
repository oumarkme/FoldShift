//! Timing of each stage (run with --ignored --nocapture). Not part of the normal test run.
use psplot_core::compare::{compare, Options};
use psplot_core::structure::load_zip;
use std::fs::File;
use std::time::Instant;

#[test]
#[ignore]
fn stages() {
    let l = |n: &str| {
        load_zip(
            File::open(format!(
                "{}/../../test/data/{n}.zip",
                env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap(),
        )
        .unwrap()
    };
    let t = Instant::now();
    let (r, c) = (l("hba1_ref"), l("hba1_trunc"));
    eprintln!("load zips: {:?}", t.elapsed());
    let t = Instant::now();
    let cmp = compare(&r, &c, &Options::default()).unwrap();
    eprintln!("compare: {:?}", t.elapsed());
    for dpi in [psplot_render::PREVIEW_DPI, psplot_render::PRINT_DPI] {
        for (ki, k) in [
            ("overlay", psplot_render::Kind::Overlay),
            ("distance", psplot_render::Kind::Distance),
        ] {
            let t = Instant::now();
            let img = psplot_render::figure(&cmp, k, &psplot_render::Style::default(), dpi);
            let t1 = t.elapsed();
            let t = Instant::now();
            let j = psplot_render::encode_jpeg(&img, dpi);
            eprintln!(
                "dpi {dpi} {}: figure {t1:?} ({}x{}), jpeg {:?} ({} KB)",
                ki,
                img.w,
                img.h,
                t.elapsed(),
                j.len() / 1024
            );
        }
    }
    let t = Instant::now();
    let p = psplot_render::profile_figure(&cmp, psplot_render::PRINT_DPI);
    eprintln!("profile figure: {:?}", t.elapsed());
    let t = Instant::now();
    psplot_render::encode_jpeg(&p, psplot_render::PRINT_DPI);
    eprintln!("profile jpeg: {:?}", t.elapsed());
}
