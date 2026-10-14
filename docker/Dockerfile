# HUSH, as a single self-contained image: the web GUI, the HTTP API and the CLI,
# with bowtie1 already inside.
#
# HUSH is distributed to users rather than hosted, so this image is the unit a
# user runs. One port serves the GUI and the API together, and the CLI is
# reachable in the same image:
#
#   docker run -p 8000:8000 -v hush_ref:/data/reference IMAGE    # GUI + API
#   docker run --rm -v "$PWD:/work" -w /work IMAGE hush --help   # CLI
#
# The transcriptome off-target reference is NOT baked in: it is multi-gigabyte
# and tracks a live RefSeq release. It is fetched once into the mounted volume
# with 'hush reference build'. See docker/README.md.

# --- Stage 1: build the web client ----------------------------------------- #
# Pinned to BUILDPLATFORM: the bundle is just files, identical on every target,
# so building it natively keeps the arm64 image off QEMU for this whole stage.
FROM --platform=$BUILDPLATFORM node:22-alpine AS web

WORKDIR /web
# Install from the lockfile first so this layer is cached until deps change.
COPY web/package.json web/package-lock.json ./
RUN npm ci
COPY web/ ./
RUN npm run build

# --- Stage 2: runtime ------------------------------------------------------ #
# bookworm is pinned deliberately: it is the suite where bowtie1 (1.3.1) is
# packaged for both amd64 and arm64, which is what makes a multi-arch image
# possible without building the aligner ourselves.
FROM python:3.12-slim-bookworm AS runtime

# bowtie1 provides `bowtie` (search) and `bowtie-build` (index). It is the only
# non-pip dependency HUSH has.
RUN apt-get update \
    && apt-get install -y --no-install-recommends bowtie \
    && rm -rf /var/lib/apt/lists/* \
    && bowtie --version | head -1 \
    && bowtie-build --version | head -1

RUN useradd --create-home --uid 10001 hush

WORKDIR /app

# Dependencies from the lockfile first, in their own layer: it is cached until
# requirements.txt changes, and it pins the environment to exactly what was
# tested rather than re-resolving ranges at build time.
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY pyproject.toml README.md LICENSE ./
COPY core/ core/
COPY msa/ msa/
COPY shrna/ shrna/
COPY service/ service/
COPY app/ app/
COPY api/ api/
COPY data/ data/
RUN pip install --no-cache-dir --no-deps .

COPY --from=web /web/dist /app/web

ENV HUSH_WEB_DIR=/app/web \
    HUSH_REFERENCE_DIR=/data/reference \
    HUSH_API_HOST=0.0.0.0 \
    HUSH_API_PORT=8000

# Pre-create the reference dir owned by the non-root run user so a fresh named
# volume inherits that ownership; otherwise the mount is root-owned and the
# index build fails with a permission error (EACCES).
RUN mkdir -p /data/reference && chown -R hush:hush /data
USER hush

# The built index and its manifest land here; mount a volume to keep them.
VOLUME ["/data/reference"]
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8000/health', timeout=3).status == 200 else 1)"

# Serve the GUI and API. Override with a command to use the CLI instead.
CMD ["python", "-m", "api"]
