use psplot_core::compare::{compare, Comparison, Options};
use psplot_core::structure::load_zip;
use psplot_render::raster::{dot, sub};
use psplot_render::scene::{camera, Camera};
use psplot_render::Style;
use std::fs::File;

fn comparison() -> Comparison {
    let load = |n: &str| {
        let path = format!("{}/../../test/data/{n}.zip", env!("CARGO_MANIFEST_DIR"));
        load_zip(File::open(path).unwrap()).unwrap()
    };
    compare(
        &load("hba1_ref"),
        &load("hba1_del_hinge"),
        &Options::default(),
    )
    .unwrap()
}

/// Picture position (pixels) of a point.
fn pixel(cam: &Camera, p: [f64; 3]) -> [f32; 2] {
    let d = sub([p[0] as f32, p[1] as f32, p[2] as f32], cam.center);
    [
        (dot(cam.axes[0], d) - cam.minxy[0]) * cam.scale,
        (-dot(cam.axes[1], d) - cam.minxy[1]) * cam.scale,
    ]
}

fn inside(cam: &Camera, p: [f32; 2]) -> bool {
    p[0] >= 0.0 && p[0] <= cam.w as f32 && p[1] >= 0.0 && p[1] <= cam.h as f32
}

#[test]
fn viewer_frame_keeps_its_size_and_holds_every_angle() {
    let c = comparison();
    let at = |yaw: f32, pitch: f32, roll: f32| {
        let st = Style {
            yaw,
            pitch,
            roll,
            steady: true,
            ..Style::default()
        };
        camera(&c, 600, &st)
    };
    let first = at(0.0, 0.0, 0.0);
    for yaw in [-135.0, -90.0, 0.0, 40.0, 90.0, 180.0] {
        for pitch in [-90.0, -30.0, 0.0, 60.0, 90.0] {
            for roll in [0.0, 75.0, -160.0] {
                let cam = at(yaw, pitch, roll);
                assert_eq!((cam.w, cam.h), (600, 600));
                assert_eq!(cam.scale, first.scale);
                assert_eq!(cam.center, first.center);
                for &p in c.ref_model.ca.iter().chain(&c.case_fit) {
                    assert!(inside(&cam, pixel(&cam, p)), "{yaw} {pitch} {roll}");
                }
            }
        }
    }
}

#[test]
fn roll_turns_the_picture_counterclockwise() {
    let c = comparison();
    let cam = |roll| {
        let st = Style {
            roll,
            steady: true,
            ..Style::default()
        };
        camera(&c, 600, &st)
    };
    let (a, b) = (cam(0.0), cam(90.0));
    // a point to the right of the middle moves above it
    for &p in &c.ref_model.ca {
        let (u, v) = (pixel(&a, p), pixel(&b, p));
        assert!((v[0] - u[1]).abs() < 0.01 && (v[1] - (600.0 - u[0])).abs() < 0.01);
    }
}

#[test]
fn zoomed_figure_shows_the_viewer_frame_without_margins() {
    let c = comparison();
    let st = Style {
        yaw: 30.0,
        zoom: 2.0,
        ..Style::default()
    };
    let figure = camera(&c, 600, &st);
    let viewer = camera(&c, 600, &Style { steady: true, ..st });
    // never smaller than on screen, and the middle of the viewer stays in the picture
    assert!(figure.scale >= viewer.scale * 0.999);
    let middle = [0, 1, 2].map(|k| f64::from(viewer.center[k]));
    assert!(inside(&figure, pixel(&figure, middle)));
}
