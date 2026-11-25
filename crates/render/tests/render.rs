use psplot_core::compare::{compare, Options};
use psplot_core::structure::load_zip;
use std::fs::File;

/// Folder for the example figures: PSPLOT_OUT, or the system temp folder.
fn out_dir() -> String {
    std::env::var("PSPLOT_OUT").unwrap_or_else(|_| std::env::temp_dir().display().to_string())
}

#[test]
#[ignore]
fn draw_examples() {
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
    let r = l("hba1_ref");
    for (tag, n) in [("del_hinge", "hba1_del_hinge"), ("trunc", "hba1_trunc")] {
        let c = compare(&r, &l(n), &Options::default()).unwrap();
        for (k, name) in [
            (psplot_render::Kind::Overlay, "overlay"),
            (psplot_render::Kind::Distance, "distance"),
            (psplot_render::Kind::Plddt, "confidence"),
        ] {
            let t = std::time::Instant::now();
            let img = psplot_render::figure(
                &c,
                k,
                &psplot_render::Style::default(),
                psplot_render::PRINT_DPI,
            );
            eprintln!("{tag} {name}: {}x{} in {:?}", img.w, img.h, t.elapsed());
            std::fs::write(
                format!("{}/rs_{tag}_{name}.jpeg", out_dir()),
                psplot_render::encode_jpeg(&img, psplot_render::PRINT_DPI),
            )
            .unwrap();
        }
        let p = psplot_render::profile_figure(&c, psplot_render::PRINT_DPI);
        std::fs::write(
            format!("{}/rs_{tag}_profile.jpeg", out_dir()),
            psplot_render::encode_jpeg(&p, psplot_render::PRINT_DPI),
        )
        .unwrap();
    }
}
