# =============================================================================
# Structure Predictor – Docker image
# Packages the full prediction pipeline (Boltz-2 + MMseqs2 + Biopython) into a
# self-contained GPU-capable container.
#
# Build:
#   docker build -t biochorl/structure-predictor:latest .
#
# Minimum host requirements:
#   • NVIDIA driver >= 525.60.13
#   • nvidia-container-toolkit
# =============================================================================

FROM nvidia/cuda:12.6.0-runtime-ubuntu22.04

LABEL maintainer="biochorl" \
      description="Protein structure retrieval and de-novo prediction pipeline" \
      version="2.0"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    # Trick: set CONDA_PREFIX so the script's conda-activation check passes
    # without modification – it believes the correct env is already active.
    CONDA_PREFIX=/home/marco/miniconda3/envs/Boltz-2

# ── System packages ──────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-dev \
        wget ca-certificates \
    && ln -sf /usr/bin/python3 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*

# ── Python packages ──────────────────────────────────────────────────────────
# Install PyTorch (CUDA 12.6) first to pin the exact build, then Boltz + deps.
RUN pip3 install --no-cache-dir \
        torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126 \
    && pip3 install --no-cache-dir \
        'boltz>=2.1' \
        biopython==1.84 \
        requests

# ── MMseqs2 (CPU-AVX2 static build) ─────────────────────────────────────────
# CPU build is sufficient for single-sequence searches and avoids GPU
# compatibility issues inside the container.
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
    && tar xzf mmseqs-linux-avx2.tar.gz \
    && cp mmseqs/bin/mmseqs /usr/local/bin/ \
    && rm -rf mmseqs mmseqs-linux-avx2.tar.gz

# ── Application ──────────────────────────────────────────────────────────────
COPY structure_predictor.py /app/structure_predictor.py

WORKDIR /workspace

ENTRYPOINT ["python3", "/app/structure_predictor.py"]
