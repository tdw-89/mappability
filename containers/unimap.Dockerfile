# UniMap.jl + samtools: the image for BUILD_LOCI, SIMULATE_AND_BATCH and
# CAT_AND_SCORE. ALIGN_BATCH uses a stock bowtie2+samtools biocontainer instead.
#
# Published as ghcr.io/tdw-89/unimap (params.container_unimap). GHCR serves
# Docker on AWS Batch and Apptainer on Unity from the same URL.
#
# Build for linux/amd64 explicitly -- Unity and Batch are both x86_64, and an
# arm64 image built on an Apple Silicon Mac pulls fine and then fails to exec:
#
#   docker buildx build --platform linux/amd64 \
#       -f containers/unimap.Dockerfile \
#       -t ghcr.io/tdw-89/unimap:0.1.0 --push .
#
# UNIMAP_REF pins the UniMap.jl commit. Bump it and the tag together. The
# commit is recorded in simulate.toml / score.toml through `git rev-parse`,
# which is why git and the .git directory stay in the image.

# Matches julia_version in UniMap.jl's Manifest.toml, so instantiate installs
# the resolved versions instead of re-resolving. trixie, not bookworm: Julia's
# bundled zstd (used to unpack the package registry) needs libzstd >= 1.5.7 and
# picks up the system copy, which on bookworm is 1.5.4.
FROM julia:1.13.2-trixie

ARG UNIMAP_REPO=https://github.com/tdw-89/UniMap.jl.git
ARG UNIMAP_REF=8db600aa3358190e93623fac3771477191338a60

# procps: Nextflow reads task metrics with ps.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        samtools \
        git \
        ca-certificates \
        procps \
        gawk \
        tar \
    && rm -rf /var/lib/apt/lists/*

# The caches are built here and loaded on whatever instance Batch picks. Without
# a multi-target build they are specialised to the build machine's CPU and get
# rejected (and silently recompiled, in every task) on a different one. This is
# the target list of the official Julia x86_64 binaries.
ENV JULIA_CPU_TARGET="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1);x86-64-v4,-rdrnd,base(1)"

RUN git clone --quiet "${UNIMAP_REPO}" /opt/UniMap.jl \
    && git -C /opt/UniMap.jl checkout --quiet "${UNIMAP_REF}"

# Precompile UniMap explicitly, not through `using UniMap` or Pkg.precompile().
# Both of those precompile every package in UniMap's [deps] -- CairoMakie and
# GLMakie included, which only `unimap view` loads -- and when one of those fails
# in a headless build, Julia quietly marks UniMap "skip precompilation" and
# loads it uncached. compilecache builds exactly what the CLI imports.
#
# One thread and serial downloads at build time only: multithreaded Julia
# segfaults in its GC safepoint under x86 emulation (an amd64 build on Apple
# Silicon). Harmless on a native builder; tasks still get -t <cpus>.
RUN export JULIA_DEPOT_PATH=/opt/julia_depot JULIA_PKG_PRECOMPILE_AUTO=0 \
           JULIA_NUM_THREADS=1,0 JULIA_PKG_CONCURRENT_DOWNLOADS=1 \
    && julia --project=/opt/UniMap.jl -e 'using Pkg; Pkg.instantiate(); \
           Base.compilecache(Base.identify_package("UniMap"))' \
    && chmod -R a+rX /opt/julia_depot /opt/UniMap.jl

# Tasks run as the calling user (docker -u, Apptainer), who cannot write to
# /opt. A writable depot first, the read-only precompiled one second: Julia
# writes logs and any new caches to the first and loads packages from either.
ENV JULIA_DEPOT_PATH=/tmp/julia_depot:/opt/julia_depot
ENV UNIMAP_PROJECT=/opt/UniMap.jl

# Fail the build, not a run, if anything is missing. Runs as an unprivileged
# user to prove the precompiled depot is usable read-only: if the shipped caches
# were rejected, Julia would quietly write new ones into /tmp/julia_depot.
RUN samtools --version | head -1 \
    && git --version \
    && su nobody -s /bin/bash -c 'cd /tmp && julia -t 1,0 --project=/opt/UniMap.jl -m UniMap --version' \
    && if [ -d /tmp/julia_depot/compiled ]; then \
           echo "precompiled caches were not reused:" >&2; ls -R /tmp/julia_depot/compiled >&2; exit 1; \
       fi \
    && rm -rf /tmp/julia_depot

CMD ["/bin/bash"]
