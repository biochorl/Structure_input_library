# =============================================================================
# Structure Predictor – Boltz-2 prediction container
#
# MINIMAL image: only Boltz-2 + PyTorch (CUDA 12.6).
# MMseqs2, Biopython, and all other pipeline logic run on the HOST.
# This container is invoked ONLY for the de-novo prediction step.
#
# Build:
#   docker build -t biochorl/structure-predictor:latest .
#
# Run (example):
#   docker run --rm --gpus all \
#       -v /path/to/workdir:/workspace \
#       -v ~/.cache/boltz_docker:/root/.cache \
#       biochorl/structure-predictor:latest \
#       predict /workspace/input.fasta \
#       --use_potentials --diffusion_samples 5 --use_msa_server
#
# Host requirements:
#   • NVIDIA driver (any version – the container ships its own CUDA 12.6)
#   • nvidia-container-toolkit
# =============================================================================

FROM nvidia/cuda:12.6.0-runtime-ubuntu22.04

LABEL maintainer="biochorl" \
      description="Boltz-2 de-novo structure prediction container" \
      version="2.0"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1

# ── System packages ──────────────────────────────────────────────────────────
# gcc is required by Triton to JIT-compile CUDA kernels at runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-dev \
        gcc \
        wget ca-certificates \
    && ln -sf /usr/bin/python3 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*

# ── Python packages ──────────────────────────────────────────────────────────
# Install PyTorch with CUDA 12.6 first, then Boltz with the [cuda] extra
# which pulls in cuequivariance-torch and cuequivariance-ops-torch-cu12
# (required for GPU-accelerated triangular multiplication kernels).
RUN pip3 install --no-cache-dir \
        torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126 \
    && pip3 install --no-cache-dir 'boltz[cuda]>=2.1'

WORKDIR /workspace

# The entrypoint IS the boltz CLI.
# Usage:  docker run ... IMAGE predict /workspace/input.fasta [boltz args]
ENTRYPOINT ["boltz"]
