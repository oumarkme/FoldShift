# Shared server

A lab can run FoldShift once in a container and let colleagues open it in a browser, with nothing installed on their PCs.

    docker run -p 8080:8080 -v foldshift-data:/data ghcr.io/<owner>/<repo>:latest     # then open http://localhost:8080/
    apptainer build foldshift.sif docker://ghcr.io/<owner>/<repo>:latest
    apptainer run --bind /path/to/data:/data foldshift.sif

The image is built from `Dockerfile` by the `release.yml` GitHub workflow (on a `v*` tag) and published to ghcr.io. Settings are environment variables:

| Variable | Default | Meaning |
|---|---|---|
| `PSPLOT_PASSWORD` | none | Optional shared password (any user name) |
| `PSPLOT_KEEP_DAYS` | 7 | Days a browser's history is kept |
| `PSPLOT_MAX_UPLOAD_MB` | 200 | Largest upload |
| `PSPLOT_MAX_JOBS` | 4 | Comparisons at the same time |
| `PSPLOT_BIND` | `0.0.0.0:8080` | Address and port |

Each browser has its own history; one browser cannot open the runs of another. The server speaks plain HTTP (no encryption), so outside a trusted network put it behind a reverse proxy with HTTPS (the page also works under a path such as `/foldshift/`, with the trailing slash). The container needs outbound internet for the transcript lookup. `/healthz` answers `ok`. `--server` on any build equals `PSPLOT_SERVER=1`.
