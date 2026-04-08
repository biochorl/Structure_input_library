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

## Installation

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

The script requires locally pre-built MMseqs2 databases for PDB and UniProtKB-TrEMBL. Navigate to the directory where you want to store them (ensure sufficient disk space) and run:

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

---

## Usage

```bash
python structure_predictor.py <input_file> [options]
```

### Command-line arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `input_file` | ✅ | — | Input file (`.fasta`, `.fa`, `.pdb`, or `.cif`) |
| `--min_identity` | ❌ | `99.0` | Minimum sequence identity (%) to accept a match |
| `--min_coverage` | ❌ | `99.0` | Minimum query coverage (%) to accept a match |
| `--pdb_db` | ❌ | *(see below)* | Absolute path to the local MMseqs2-formatted PDB database |
| `--uniprot_db` | ❌ | *(see below)* | Absolute path to the local MMseqs2-formatted UniProtKB-TrEMBL database |

### Examples

**Basic run with a FASTA file:**
```bash
python structure_predictor.py my_protein.fasta
```

**Specifying custom database paths (required on different machines):**
```bash
python structure_predictor.py my_protein.fasta \
    --pdb_db /data/mmseqs_dbs/PDB \
    --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL
```

**Custom identity and coverage thresholds:**
```bash
python structure_predictor.py my_protein.fasta --min_identity 95.0 --min_coverage 90.0
```

The script will automatically handle Conda environment activation and proceed with the search or prediction.

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

The PDB and UniProtKB-TrEMBL databases are searched **locally** via pre-built MMseqs2 databases (see [Installation](#installation)). AlphaFold DB is queried **remotely** via its REST API to retrieve pre-computed models.

---

## Third-party software & licenses

This pipeline integrates or calls the following tools. They are **not bundled** with this repository — users must install them separately.

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
