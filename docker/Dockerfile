# HUSH container image: Python 3.11 + bowtie1, with HUSH installed on top.
#
# One image serves both off-target roles so bowtie always lives on our side,
# never on the user's machine:
#   1. building the RefSeq off-target index (tools/build_offtarget_db.py), and
#   2. the future hosted design backend, which shells out to bowtie at run time.
#
# Because it is a Linux container, it also removes the host-OS friction of
# installing bowtie1 via bioconda on native Windows: the build runs the same way
# on any host with Docker. See docker/README.md for the build/run commands.

FROM mambaorg/micromamba:1.5.8

# Create the environment (Python + bowtie1) from the pinned spec first, so this
# heavy layer is cached and only re-run when the spec changes.
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate the base env for every following RUN, and at container start.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
RUN pip install --no-cache-dir .

# Pre-create the reference dir owned by the non-root run user so a fresh named
# volume inherits that ownership; otherwise the mount is root-owned and the
# index build fails with a permission error (EACCES).
USER root
RUN mkdir -p /data/reference && chown -R $MAMBA_USER:$MAMBA_USER /data
USER $MAMBA_USER

# The built index and its manifest land here; mount a volume to persist them
# across containers and to hand the index to the runtime backend.
VOLUME ["/data/reference"]

# Default to the CLI help. Override the command to build the index, e.g.
#   docker run --rm -v hush_ref:/data/reference hush \
#     python tools/build_offtarget_db.py --outdir /data/reference
CMD ["hush", "--help"]
