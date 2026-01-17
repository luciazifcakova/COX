FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# ---- Minimal OS deps ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl \
    pigz \
    libgomp1 \
    bzip2 \
  && rm -rf /var/lib/apt/lists/*

# ---- micromamba install ----
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# Use bash -lc so micromamba commands behave consistently
SHELL ["/bin/bash", "-lc"]

# ---- Create env for pipeline ----
RUN micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config set channel_priority strict && \
    micromamba create -y -n coi_pipeline \
      python=3.11 pip \
      pandas \
      numpy \
      cutadapt \
      nanofilt \
      nanoplot \
      multiqc \
      minimap2 \
      samtools \
      racon \
      spoa \
      medaka=2.0.1 \
    && micromamba clean -a -y
# Put the env first on PATH
ENV PATH="/opt/micromamba/envs/coi_pipeline/bin:${PATH}"

# ---- NGSpeciesID ----
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir NGSpeciesID

# ---- BLAST+ install ----
# (contains blastn, makeblastdb, blastdbcmd, etc.)
WORKDIR /tmp
RUN curl -Lso ncbi-blast-2.16.0+-x64-linux.tar.gz \
      https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar xzf ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.16.0+/bin/* /usr/local/bin/ && \
    rm -rf /tmp/ncbi-blast-2.16.0* /tmp/ncbi-blast-2.16.0+-x64-linux.tar.gz

# ---- Your pipeline scripts ----
WORKDIR /app
COPY scripts /app/scripts
RUN find /app/scripts -type f \( -name "*.sh" -o -name "*.py" \) -exec chmod +x {} \;

ENV PATH="/app/scripts:${PATH}"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN mkdir -p /data/input /data/output /data/databases

# ---- Quick sanity checks (only what you actually use) ----
RUN NanoPlot --version >/dev/null && \
    multiqc --version >/dev/null && \
    NGSpeciesID --help >/dev/null && \
    medaka_consensus -h >/dev/null && \
    minimap2 --version >/dev/null && \
    samtools --version >/dev/null && \
    blastn -version >/dev/null && \
    python -c "import numpy,pandas; print('python env ok')"
WORKDIR /data
CMD ["/bin/bash"]
