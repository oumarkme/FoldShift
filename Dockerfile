# FoldShift, server mode in one small image.
#   docker build -t foldshift .
#   docker run -p 8080:8080 -v psplot-data:/data foldshift
# Then open http://localhost:8080/ . Settings (environment variables): PSPLOT_PASSWORD (optional shared password),
# PSPLOT_KEEP_DAYS (default 7), PSPLOT_MAX_UPLOAD_MB (default 200), PSPLOT_MAX_JOBS (default 4).
FROM rust:bookworm AS build
WORKDIR /src
COPY rust-toolchain.toml Cargo.toml Cargo.lock ./
COPY crates ./crates
RUN cargo build --release --locked -p psplot

FROM debian:bookworm-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates wget \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --system --uid 10001 --create-home --home-dir /data psplot
COPY --from=build /src/target/release/FoldShift /usr/local/bin/FoldShift
ENV PSPLOT_SERVER=1 PSPLOT_HOME=/data PSPLOT_BIND=0.0.0.0:8080
USER psplot
VOLUME /data
EXPOSE 8080
HEALTHCHECK --interval=30s --timeout=3s CMD wget -qO- http://127.0.0.1:8080/healthz >/dev/null || exit 1
ENTRYPOINT ["FoldShift", "--server"]
