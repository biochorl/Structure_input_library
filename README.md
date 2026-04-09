# 🧬 Structure Predictor

Automated protein structure retrieval and *de novo* prediction pipeline.

Given a protein sequence (FASTA) or an existing structure file (PDB/CIF), the script follows a hierarchical strategy to obtain the best available 3D structure as efficiently as possible.

---

## Overview

```
Input FASTA
    │
    ├─► MMseqs2 search against local PDB database
    │       └─► Match found → download from RCSB, align & truncate → ✅ done
    │
    ├─► MMseqs2 search against local UniProtKB-TrEMBL database
    │       └─► Match found → fetch AlphaFold DB model (pLDDT ≥ 75), align & truncate → ✅ done
    │
    └─► De novo prediction with Boltz-2
            └─► Best model (confidence ≥ 0.8) → ✅ done
```

### Step-by-step

1. **Input analysis** – Determines whether the input is a sequence (FASTA) or a structure (PDB/CIF).
   - **Structure file**: extracts the first protein chain and saves it as a new PDB.
   - **FASTA file**: starts the hierarchical search described below.

2. **Hierarchical structure search** (for FASTA input):
   1. **PDB search** – Runs a local MMseqs2 search against a pre-built PDB database to find an experimentally resolved structure matching the query above the identity/coverage thresholds. If found, the structure is downloaded from RCSB PDB.
   2. **AlphaFold DB lookup** – If no PDB match is found, searches a local UniProtKB-TrEMBL database. The resulting UniProt accession is used to query the AlphaFold DB API. If a high-quality pre-computed model is available (**pLDDT ≥ 75**), it is downloaded.
   3. **De novo prediction with Boltz-2** – As a last resort, the script runs Boltz-2 locally to generate 5 candidate models and selects the one with the highest confidence score (threshold: **≥ 0.8**).

3. **Alignment & truncation** – Every candidate structure is locally re-aligned to the query sequence with Biopython. If identity and coverage thresholds are met, the model is truncated to the aligned region and saved.

---

## Two installation modes

| | **v1 – Conda** (`structure_predictor.py`) | **v2 – Docker** (`structure_predictor_docker.py`) |
|---|---|---|
| **Host requirements** | Conda, CUDA toolkit, NVIDIA driver matching your PyTorch build | Docker, nvidia-container-toolkit, NVIDIA driver ≥ 525.60.13 |
| **Best for** | Systems with full CUDA stack (Ubuntu 22.04+, recent drivers) | Older systems (e.g. Ubuntu 20.04) or when you want a reproducible, self-contained environment |
| **GPU access** | Direct (via Conda environment) | Via NVIDIA Container Toolkit |

---

## Input & Output

### Supported input files

| Format | Extensions |
|--------|-----------|
| FASTA  | `.fasta`, `.fa` |
| PDB    | `.pdb` |
| mmCIF  | `.cif` |

### Generated output files

For an input named `my_protein.fasta`:

| File / Directory | Description |
|-----------------|-------------|
| `my_protein.pdb` | Final 3D structure (single chain) |
| `my_protein_report.log` | Timestamped execution log |
| `my_protein_temp_files/` | Temporary downloads (auto-deleted) |
| `my_protein_boltz_results/` | Boltz-2 prediction output (if triggered) |

---

## Installation – v1 (Conda)

### 1. Create a Conda environment

```bash
conda create -n boltz python=3.10
conda activate boltz
```

### 2. Install dependencies

**Boltz-2** ([MIT License](https://github.com/jwohlwend/boltz/blob/main/LICENSE))

```bash
pip install 'boltz[cuda]' -U
```

**MMseqs2 – GPU build** ([GPLv3 License](https://github.com/soedinglab/MMseqs2/blob/master/LICENSE.md))

Download the pre-compiled GPU-accelerated build. This version requires glibc ≥ 2.29 and NVIDIA driver ≥ 525.60.13.

```bash
# Download and extract the MMseqs2-GPU archive
wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz
tar xvfz mmseqs-linux-gpu.tar.gz
```

After extraction, add the binary to your PATH permanently by appending it to `~/.bashrc` and reloading:

```bash
echo "export PATH=$(pwd)/mmseqs/bin/:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

**Download the MMseqs2 databases**

The script requires locally pre-built MMseqs2 databases for PDB and UniProtKB-TrEMBL. The combined databases require approximately **~120 GB** of disk space. Navigate to the directory where you want to store them and run:

```bash
mkdir -p tmp
mmseqs databases PDB PDB tmp
mmseqs databases UniProtKB/TrEMBL UniProtKB-TrEMBL tmp
# The temporary directory can be removed afterwards
rm -rf tmp
```

> **Note:** Remember the absolute path to this directory — you will need to pass it to the script via `--pdb_db` and `--uniprot_db` if the databases are not stored in the default location.

**Biopython** ([Biopython License / BSD-3-Clause](https://github.com/biopython/biopython/blob/master/LICENSE.rst))

```bash
pip install biopython
```

### 3. Script configuration

Set the `CONDA_ENV_PATH` variable at the top of `structure_predictor.py` to the absolute path of the Conda environment you just created. The script will automatically re-launch itself inside the correct environment if it is not already active.

### Usage (v1)

```bash
python structure_predictor.py <input_file> [options]
```

#### Command-line arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `input_file` | ✅ | — | Input file (`.fasta`, `.fa`, `.pdb`, or `.cif`) |
| `--min_identity` | ❌ | `99.0` | Minimum sequence identity (%) to accept a match |
| `--min_coverage` | ❌ | `99.0` | Minimum query coverage (%) to accept a match |
| `--pdb_db` | ❌ | *(see below)* | Absolute path to the local MMseqs2-formatted PDB database |
| `--uniprot_db` | ❌ | *(see below)* | Absolute path to the local MMseqs2-formatted UniProtKB-TrEMBL database |

#### Examples (v1)

```bash
# Basic run
python structure_predictor.py my_protein.fasta

# Custom database paths
python structure_predictor.py my_protein.fasta \
    --pdb_db /data/mmseqs_dbs/PDB \
    --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL

# Custom thresholds
python structure_predictor.py my_protein.fasta --min_identity 95.0 --min_coverage 90.0
```

---

## Installation – v2 (Docker)

The Docker version packages all dependencies (Boltz-2, PyTorch, MMseqs2, Biopython) into a single container image. **No Conda environment, CUDA toolkit, or Python scientific packages are needed on the host.** This is the recommended method for systems where installing the latest CUDA stack is difficult or impossible (e.g. Ubuntu 20.04 workstations with older NVIDIA drivers).

### 1. Install Docker

Follow the official instructions for your distribution: https://docs.docker.com/engine/install/

For Ubuntu 20.04:

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates curl gnupg
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu focal stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io
```

Add your user to the `docker` group so you can run containers without `sudo`:

```bash
sudo usermod -aG docker $USER
newgrp docker    # or log out and back in
```

### 2. Install NVIDIA driver (≥ 525.60.13)

The Docker container runs its own CUDA toolkit internally, but the host still needs a compatible NVIDIA GPU driver. The minimum required version is **525.60.13**.

For Ubuntu 20.04, install a supported driver from the NVIDIA PPA:

```bash
sudo add-apt-repository -y ppa:graphics-drivers/ppa
sudo apt-get update
# Install driver 535 (or any version >= 525)
sudo apt-get install -y nvidia-driver-535
sudo reboot
```

After rebooting, verify with:

```bash
nvidia-smi
```

### 3. Install NVIDIA Container Toolkit

This enables Docker containers to access the host GPU:

```bash
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey \
  | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list \
  | sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' \
  | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list > /dev/null
sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

Verify GPU access inside Docker:

```bash
docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu22.04 nvidia-smi
```

### 4. Download the MMseqs2 databases

Same as the Conda version — see above. The databases live on the host and are bind-mounted into the container at runtime.

### 5. Build the Docker image

```bash
cd structure-predictor/
docker build -t biochorl/structure-predictor:latest .
```

Or let the wrapper script build it for you with the `--build` flag (see below).

### Usage (v2)

```bash
python3 structure_predictor_docker.py <input_file> [options]
```

> **Note:** The wrapper script only requires Python 3 (any version ≥ 3.6) and Docker on the host. No scientific packages are needed.

#### Command-line arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `input_file` | ✅ | — | Input file (`.fasta`, `.fa`, `.pdb`, or `.cif`) |
| `--min_identity` | ❌ | `99.0` | Minimum sequence identity (%) |
| `--min_coverage` | ❌ | `99.0` | Minimum query coverage (%) |
| `--pdb_db` | ✅ | — | Absolute path to the local MMseqs2-formatted PDB database |
| `--uniprot_db` | ✅ | — | Absolute path to the local MMseqs2-formatted UniProtKB-TrEMBL database |
| `--docker_image` | ❌ | `biochorl/structure-predictor:latest` | Docker image name |
| `--boltz_cache` | ❌ | `~/.cache/boltz_docker` | Host path for persistent Boltz model weights |
| `--build` | ❌ | — | Build/rebuild the Docker image before running |

#### Examples (v2)

```bash
# Build the image and run in one step
python3 structure_predictor_docker.py my_protein.fasta --build \
    --pdb_db /data/mmseqs_dbs/PDB \
    --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL

# Run (image already built)
python3 structure_predictor_docker.py my_protein.fasta \
    --pdb_db /data/mmseqs_dbs/PDB \
    --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL

# Custom thresholds and persistent cache location
python3 structure_predictor_docker.py my_protein.fasta \
    --min_identity 95.0 --min_coverage 90.0 \
    --pdb_db /data/mmseqs_dbs/PDB \
    --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL \
    --boltz_cache /data/boltz_weights
```

#### How it works

The wrapper script (`structure_predictor_docker.py`) translates your host paths into Docker volume mounts and runs the original `structure_predictor.py` inside the container:

```
Host                                    Container
─────────────────────────────────────   ──────────────────────
./my_protein.fasta                  →   /workspace/my_protein.fasta
/data/mmseqs_dbs/  (PDB files)      →   /db/pdb/    (read-only)
/data/mmseqs_dbs/  (UniProt files)  →   /db/uniprot/ (read-only)
~/.cache/boltz_docker               →   /root/.cache (read-write)
```

Output files (`my_protein.pdb`, `my_protein_report.log`, etc.) are written to the same host directory as the input file.

---

## Example files

The `examples/` directory contains sample input files you can use to test the pipeline:

| File | Description |
|------|-------------|
| `Input_PDB.fasta` | Sequence with a known PDB match |
| `Input_AlphaFoldDB.fasta` | Sequence with a match on AlphaFold DB |
| `Input_Boltz_prediction.fasta` | Sequence requiring *de novo* prediction |

---

## Databases

The pipeline queries the following public biological databases:

| Database | URL | Description |
|----------|-----|-------------|
| [RCSB Protein Data Bank (PDB)](https://www.rcsb.org/) | https://www.rcsb.org/ | Archive of experimentally determined 3D structures of proteins, nucleic acids, and complex assemblies |
| [UniProt / UniProtKB-TrEMBL](https://www.uniprot.org/) | https://www.uniprot.org/ | Comprehensive, high-quality protein sequence and functional annotation resource |
| [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) | https://alphafold.ebi.ac.uk/ | Database of protein structure predictions by AlphaFold, covering the UniProt proteomes |

The PDB and UniProtKB-TrEMBL databases are searched **locally** via pre-built MMseqs2 databases (see [Installation](#installation--v1-conda)). AlphaFold DB is queried **remotely** via its REST API to retrieve pre-computed models.

---

## Third-party software & licenses

This pipeline integrates or calls the following tools. They are **not bundled** with this repository — users must install them separately (or use the Docker image which includes them all).

| Software | License | Role in the pipeline |
|----------|---------|---------------------|
| [Boltz-2](https://github.com/jwohlwend/boltz) | MIT | *De novo* structure prediction |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | GPLv3 | Local sequence database search |
| [Biopython](https://github.com/biopython/biopython) | Biopython License / BSD-3-Clause | Sequence I/O, alignment, PDB parsing |

This project's own code is released under the [MIT License](LICENSE).

> **License compatibility note:** MMseqs2 is licensed under GPLv3. This project calls MMseqs2 as an external command-line tool (not linked as a library), so the GPLv3 license of MMseqs2 does not propagate to this project's MIT-licensed code. Users must comply with the GPLv3 when using or distributing MMseqs2 itself.

---

## Citation

If you use this tool in your research, please cite the underlying software and databases:

```bibtex
@article{wohlwend2024boltz1,
  title   = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  author  = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Reveiz, Mateo
             and Leidal, Ken and Swanson, Wojtek and Gilmer, Justin and Ronning, Bowen
             and Husic, Brooke and Durairaj, Jiarui and others},
  journal = {bioRxiv},
  year    = {2024}
}

@article{steinegger2017mmseqs2,
  title   = {MMseqs2 enables sensitive protein sequence searching for the
             analysis of massive data sets},
  author  = {Steinegger, Martin and S{\"o}ding, Johannes},
  journal = {Nature Biotechnology},
  volume  = {35},
  pages   = {1026--1028},
  year    = {2017}
}

@article{berman2000pdb,
  title   = {The Protein Data Bank},
  author  = {Berman, Helen M and Westbrook, John and Feng, Zukang and Gilliland, Gary
             and Bhat, T N and Weissig, Helge and Shindyalov, Ilya N and Bourne, Philip E},
  journal = {Nucleic Acids Research},
  volume  = {28},
  number  = {1},
  pages   = {235--242},
  year    = {2000}
}

@article{uniprot2023,
  title   = {{UniProt}: the Universal Protein Knowledgebase in 2023},
  author  = {{The UniProt Consortium}},
  journal = {Nucleic Acids Research},
  volume  = {51},
  number  = {D1},
  pages   = {D523--D531},
  year    = {2023}
}

@article{varadi2022alphafolddb,
  title   = {{AlphaFold Protein Structure Database}: massively expanding the
             structural coverage of protein-sequence space with high-accuracy models},
  author  = {Varadi, Mihaly and Anyango, Stephen and Deshpande, Mandar and Nair, Sreenath
             and Natassia, Cindy and Yordanova, Galabina and Yuan, David and Stroe, Oana
             and Wood, Gemma and Laydon, Agata and others},
  journal = {Nucleic Acids Research},
  volume  = {50},
  number  = {D1},
  pages   = {D439--D444},
  year    = {2022}
}
```
