#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structure Predictor – Docker edition (v2)

Full prediction pipeline that runs natively on the host.  Only the
Boltz-2 de-novo prediction step is offloaded to a Docker container
with GPU support.

Host requirements:
    • Python 3.8+  (with biopython and requests)
    • MMseqs2      (CPU build is fine)
    • Docker       (≥ 19.03)
    • nvidia-container-toolkit  (for GPU access inside the container)
    • NVIDIA driver              (any version supported by your GPU)
"""

import os
import sys
import argparse
import subprocess
import shutil
import requests
import json
import datetime
import gzip
from glob import glob

try:
    from Bio import SeqIO, Align
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.MMCIFParser import MMCIFParser
except ImportError:
    print("ERROR: Biopython is not installed. Install with 'pip install biopython'",
          file=sys.stderr)
    sys.exit(1)

DOCKER_IMAGE = "biochorl/structure-predictor:latest"
RESIDUE_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V"
}

# ── Logging ──────────────────────────────────────────────────────────────────

class Logger:
    def __init__(self, log_file_path):
        self.terminal = sys.stdout
        self.log_file = open(log_file_path, 'w')

    def write(self, message):
        self.terminal.write(message)
        if message.strip():
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            self.log_file.write(f"[{timestamp}] {message.rstrip()}\n")

    def flush(self):
        self.terminal.flush()
        self.log_file.flush()

    def __del__(self):
        if not self.log_file.closed:
            self.log_file.close()
        sys.stdout = self.terminal

# ── PDB selection helpers ────────────────────────────────────────────────────

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

class ResidueSelect(Select):
    def __init__(self, chain_id, start, end):
        self.chain_id = chain_id
        self.start = start
        self.end = end

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue):
        if residue.get_id()[0] != ' ': return False
        return self.start <= residue.get_id()[1] <= self.end

# ── Core functions (run on the HOST) ─────────────────────────────────────────

def extract_seq_from_pdb(pdb_path, chain_id):
    print(f"🧬 Extracting sequence from chain '{chain_id}'…")
    if pdb_path.endswith(('.cif', '.mmcif')):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    try:
        try: structure = parser.get_structure("s", pdb_path)
        except (UnicodeDecodeError, IsADirectoryError):
            with gzip.open(pdb_path, 'rt') as f:
                if pdb_path.endswith(('.cif', '.mmcif')):
                    gz_parser = MMCIFParser(QUIET=True)
                else:
                    gz_parser = PDBParser(QUIET=True)
                structure = gz_parser.get_structure("s", f)

        chain = structure[0][chain_id]
        seq = [RESIDUE_MAP.get(r.get_resname().strip().upper(), 'X')
               for r in chain.get_residues() if r.get_id()[0] == ' ']
        return "".join(seq)
    except Exception as e:
        print(f"ERROR extracting PDB sequence: {e}")
        return None

def align_and_verify(query_seq, target_seq, min_identity, min_coverage):
    print("\n" + "="*35)
    print("🔬  STARTING PAIRWISE ALIGNMENT  🔬")
    print("="*35)
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1

    alignments = aligner.align(query_seq, target_seq)
    if not alignments: return False, 0, 0, None, None

    best_alignment = alignments[0]
    print("\n[Best alignment]:\n" + str(best_alignment))

    aligned_q, aligned_t = best_alignment
    matches, alignment_len_no_gaps, query_aligned_len = 0, 0, 0
    t_start_index, t_end_index, t_res_count = -1, -1, 0

    for i, (q_char, t_char) in enumerate(zip(aligned_q, aligned_t)):
        if q_char != '-': query_aligned_len += 1
        if t_char != '-':
            if t_start_index == -1: t_start_index = t_res_count
            t_end_index = t_res_count
            t_res_count += 1
        if q_char != '-' and t_char != '-':
            alignment_len_no_gaps += 1
            if q_char == t_char: matches += 1

    identity = (matches / alignment_len_no_gaps) * 100 if alignment_len_no_gaps > 0 else 0
    coverage = (query_aligned_len / len(query_seq)) * 100

    t_start_pdb, t_end_pdb = (t_start_index + 1, t_end_index + 1) if t_start_index != -1 else (None, None)

    print("\n" + "-"*35)
    print("📊 Alignment results:")
    print(f"   - Identity: {identity:.2f}%")
    print(f"   - Coverage: {coverage:.2f}%")
    print(f"   - Aligned PDB region: residues {t_start_pdb} – {t_end_pdb}")
    print("-"*35 + "\n")

    is_match = identity >= min_identity and coverage >= min_coverage
    if is_match: return True, identity, coverage, t_start_pdb, t_end_pdb
    return False, identity, coverage, None, None

def truncate_and_save_pdb(input_pdb_path, output_pdb_path, chain_id, start, end):
    print(f"✂️  Truncating model (residues {start}–{end}, chain {chain_id})…")
    if input_pdb_path.endswith(('.cif', '.mmcif')):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    try:
        try: structure = parser.get_structure("full_model", input_pdb_path)
        except (UnicodeDecodeError, IsADirectoryError):
            with gzip.open(input_pdb_path, 'rt') as f:
                if input_pdb_path.endswith(('.cif', '.mmcif')):
                    gz_parser = MMCIFParser(QUIET=True)
                else:
                    gz_parser = PDBParser(QUIET=True)
                structure = gz_parser.get_structure("full_model", f)
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_path, ResidueSelect(chain_id, start, end))
        print(f"✔️ Truncated model saved: {output_pdb_path}")
        return True
    except Exception as e:
        print(f"ERROR during truncation: {e}", file=sys.stderr)
        shutil.copy(input_pdb_path, output_pdb_path)
        return False

def process_pdb_file(pdb_path, output_path):
    print(f"📖 Processing structure file: {pdb_path}")
    if pdb_path.endswith(('.cif', '.mmcif')):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("s", pdb_path)
        chain = next(structure[0].get_chains())
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_path, ChainSelect(chain.id))
        print(f"✔️ Chain '{chain.id}' saved to: {output_path}")
        return True
    except Exception as e:
        print(f"ERROR processing PDB/CIF: {e}", file=sys.stderr)
        return False

def check_mmseqs_availability():
    if not shutil.which("mmseqs"):
        sys.exit("ERROR: 'mmseqs' not found in PATH. Install MMseqs2.")
    print("✔️ mmseqs found.")

def parse_fasta(fasta_path):
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq).upper(), len(record.seq)
    except Exception as e:
        sys.exit(f"ERROR reading FASTA: {e}")

def run_mmseqs_search(temp_fasta, query_len, db_path, temp_dir,
                      min_identity, min_coverage, extra_params=None):
    db_name = os.path.basename(db_path)
    print(f"\n🔎 Running MMseqs2 search against '{db_name}'…")
    out_m8 = os.path.join(temp_dir, f"search_{db_name}.m8")
    tmp_mmseqs = os.path.join(temp_dir, f"mmseqs_tmp_{db_name}")
    os.makedirs(tmp_mmseqs, exist_ok=True)
    command = ["mmseqs", "easy-search", temp_fasta, db_path, out_m8,
               tmp_mmseqs, "--format-output",
               "target,pident,alnlen,qlen,tlen,tstart,tend"]
    if extra_params:
        command.extend(extra_params)
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        if not os.path.exists(out_m8):
            print("MMseqs2 produced no results."); return None
        with open(out_m8, 'r') as f:
            hits = f.read().strip().split('\n')
        if not hits or not hits[0]:
            print("MMseqs2 produced no results."); return None
        for line in hits:
            parts = line.split('\t')
            if len(parts) < 7: continue
            sseqid, pident, align_len, qlen, slen, sstart, send = parts
            coverage = (int(align_len) / query_len) * 100
            if float(pident) >= min_identity and coverage >= min_coverage:
                print(f"✔️ Match found in {db_name}: {sseqid.strip()}")
                return (sseqid.strip(), int(sstart), int(send))
        print(f"\n❌ No match in {db_name} passed the thresholds.")
        return None
    except Exception as e:
        print(f"ERROR MMseqs2: {e}", file=sys.stderr)
        return None

def download_pdb_file(pdb_id, output_path, source="rcsb"):
    url = {"rcsb": f"https://files.rcsb.org/download/{pdb_id}.pdb",
           "alphafold": f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb"
           }.get(source)
    print(f"⬇️  Downloading PDB from {url}…")
    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(output_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
            print(f"✔️ Saved: {output_path}")
            return True
    except requests.exceptions.RequestException as e:
        print(f"ERROR: Download failed. {e}", file=sys.stderr)
        return False

def get_alphafold_structure(uniprot_id, output_path):
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"\n🔎 Querying AlphaFold DB API: {api_url}…")
    try:
        res = requests.get(api_url, timeout=60)
        if res.status_code != 200:
            print(f"❌ API request failed (status {res.status_code})",
                  file=sys.stderr)
            return False
        try:
            prediction_data = res.json()
        except requests.exceptions.JSONDecodeError:
            print("❌ API response is not valid JSON.", file=sys.stderr)
            return False
        if not prediction_data:
            print("❌ API returned empty response.", file=sys.stderr)
            return False
        entry_data = prediction_data[0]
        model_url = entry_data.get('pdbUrl') or entry_data.get('cifUrl')
        avg_plddt = entry_data.get('globalMetricValue')
        if not model_url or avg_plddt is None:
            print("❌ No model URL or pLDDT score in response.",
                  file=sys.stderr)
            return False
        print(f"📊 AlphaFold DB quality check: pLDDT={avg_plddt:.2f}")
        if avg_plddt >= 75:
            print("✔️ Model meets quality criteria (pLDDT ≥ 75).")
            print(f"⬇️  Downloading model from {model_url}…")
            with requests.get(model_url, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(output_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)
                print(f"✔️ Saved: {output_path}")
                return True
        else:
            print("❌ Model does not meet quality threshold.")
            return False
    except Exception as e:
        print(f"ERROR querying AlphaFold DB: {e}", file=sys.stderr)
        return False

# ── Boltz prediction via Docker ──────────────────────────────────────────────

def run_boltz(fasta_path, final_output_path, boltz_output_dir,
              docker_image, boltz_cache):
    """Run Boltz-2 de-novo prediction *inside a Docker container*.

    Only this function uses Docker.  Everything else in the pipeline
    runs natively on the host.
    """
    print(f"\n🚀 Starting Boltz-2 prediction (via Docker)…")
    if os.path.exists(boltz_output_dir):
        shutil.rmtree(boltz_output_dir)
    os.makedirs(boltz_output_dir, exist_ok=True)
    os.makedirs(boltz_cache, exist_ok=True)

    # Prepare input FASTA with Boltz-compatible header
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
    except Exception as e:
        print(f"ERROR reading FASTA for Boltz: {e}", file=sys.stderr)
        return False

    temp_fasta_name = "temp_input.fasta"
    temp_fasta_host = os.path.join(boltz_output_dir, temp_fasta_name)
    with open(temp_fasta_host, "w") as f:
        f.write(f">A|protein|\n{record.seq}\n")

    # ── Build docker run command ─────────────────────────────────────────
    abs_output_dir = os.path.abspath(boltz_output_dir)
    abs_cache      = os.path.abspath(boltz_cache)

    docker_cmd = [
        "docker", "run", "--rm",
        "--gpus", "all",
        "-v", f"{abs_output_dir}:/workspace",
        "-v", f"{abs_cache}:/root/.cache",
        "-w", "/workspace",
        docker_image,
        # ── boltz CLI arguments ──
        "predict",
        f"/workspace/{temp_fasta_name}",
        "--use_potentials",
        "--diffusion_samples", "5",
        "--use_msa_server",
    ]

    print(f"🐳 Docker command:\n   {' '.join(docker_cmd)}\n")

    try:
        result = subprocess.run(docker_cmd, timeout=3600)
        if result.returncode != 0:
            print("❌ Boltz Docker container exited with non-zero status.",
                  file=sys.stderr)
            return False
    except subprocess.TimeoutExpired:
        print("❌ Boltz prediction timed out (1 h).", file=sys.stderr)
        return False
    except FileNotFoundError:
        print("❌ 'docker' command not found. Is Docker installed?",
              file=sys.stderr)
        return False

    # ── Parse results (files are on the HOST via the bind mount) ─────────
    fasta_stem = os.path.splitext(temp_fasta_name)[0]   # "temp_input"
    prediction_base_dir = os.path.join(
        boltz_output_dir,
        f"boltz_results_{fasta_stem}", "predictions", fasta_stem)

    score_files = glob(os.path.join(prediction_base_dir, "confidence_*.json"))
    if not score_files:
        print("❌ No confidence score files found.", file=sys.stderr)
        return False

    ranked_models = sorted(
        [(json.load(open(f)).get('confidence_score'),
          f.replace("confidence_", "").replace(".json", ".cif"))
         for f in score_files],
        reverse=True)

    best_score, best_model_file = ranked_models[0]
    print(f"\n📊 Best confidence score: {best_score:.4f}")
    if best_score >= 0.8:
        print("✔️ Model meets quality criteria.")
        return process_pdb_file(best_model_file, final_output_path)
    else:
        print("❌ INSUFFICIENT QUALITY.")
        return False

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Structure Predictor – Docker edition (v2). "
                    "Searches run on the host; Boltz-2 prediction runs "
                    "inside a Docker container with GPU support.")
    parser.add_argument("input_file",
                        help="Input file (.fasta, .fa, .pdb, or .cif).")
    parser.add_argument("--min_identity", type=float, default=99.0,
                        help="Minimum sequence identity (%%) [default: 99.0].")
    parser.add_argument("--min_coverage", type=float, default=99.0,
                        help="Minimum query coverage (%%) [default: 99.0].")
    parser.add_argument("--pdb_db", type=str,
                        default="/media/marco/Database/Test_Database/PDB",
                        help="Path to local MMseqs2-formatted PDB database.")
    parser.add_argument("--uniprot_db", type=str,
                        default="/media/marco/Database/Test_Database/UniProtKB-TrEMBL",
                        help="Path to local MMseqs2-formatted UniProtKB "
                             "database.")
    parser.add_argument("--docker_image", type=str, default=DOCKER_IMAGE,
                        help=f"Docker image for Boltz-2 "
                             f"[default: {DOCKER_IMAGE}].")
    parser.add_argument("--boltz_cache", type=str,
                        default=os.path.expanduser("~/.cache/boltz_docker"),
                        help="Host directory for persistent Boltz model "
                             "weights [default: ~/.cache/boltz_docker].")
    args = parser.parse_args()

    input_path = os.path.abspath(args.input_file)
    base_name  = os.path.splitext(os.path.basename(input_path))[0]
    output_dir = os.path.dirname(input_path)
    output_pdb = os.path.join(output_dir, f"{base_name}.pdb")
    log_file   = os.path.join(output_dir, f"{base_name}_report.log")

    sys.stdout = Logger(log_file)

    print("✔️ Structure Predictor – Docker edition (v2)")

    temp_dir = os.path.join(output_dir, f"{base_name}_temp_files")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)

    print(f"Processing: {os.path.basename(input_path)}")

    # ── PDB / CIF input ─────────────────────────────────────────────────
    if input_path.endswith(('.pdb', '.cif')):
        process_pdb_file(input_path,
                         os.path.join(output_dir, f"{base_name}_chain1.pdb"))

    # ── FASTA input ──────────────────────────────────────────────────────
    elif input_path.endswith(('.fasta', '.fa')):
        check_mmseqs_availability()
        sequence, seq_len = parse_fasta(input_path)
        print(f"Sequence read: {seq_len} aa.")

        temp_query = os.path.join(temp_dir, "query.fasta")
        with open(temp_query, "w") as f:
            f.write(f">query\n{sequence}\n")

        pdb_db     = os.path.abspath(args.pdb_db)
        uniprot_db = os.path.abspath(args.uniprot_db)

        # Step 1 – PDB search (HOST)
        pdb_hit = run_mmseqs_search(temp_query, seq_len, pdb_db,
                                    temp_dir, args.min_identity,
                                    args.min_coverage)
        if pdb_hit:
            sseqid, _, _ = pdb_hit
            parts = sseqid.replace('_', '|').split('|')
            if len(parts) >= 2:
                pdb_code = parts[-2][:4].lower()
                chain_id = parts[-1]
                if len(pdb_code) == 4:
                    tmp_pdb = os.path.join(temp_dir, f"{pdb_code}.pdb")
                    if download_pdb_file(pdb_code, tmp_pdb, source="rcsb"):
                        pdb_seq = extract_seq_from_pdb(tmp_pdb, chain_id)
                        if pdb_seq:
                            ok, *_ , start, end = align_and_verify(
                                sequence, pdb_seq,
                                args.min_identity, args.min_coverage)
                            if ok:
                                truncate_and_save_pdb(tmp_pdb, output_pdb,
                                                     chain_id, start, end)
                                print("\n🎉 Done!")
                                shutil.rmtree(temp_dir)
                                sys.exit(0)

        # Step 2 – AlphaFold DB lookup (HOST)
        # Use faster search parameters for the large UniProt database:
        # -s 4 reduces sensitivity (default 5.7) for ~3x speedup,
        # --max-seqs 100 limits prefilter results since we only need the top hit.
        uniprot_fast_params = ["-s", "4", "--max-seqs", "100"]
        uniprot_hit = run_mmseqs_search(temp_query, seq_len, uniprot_db,
                                        temp_dir, args.min_identity,
                                        args.min_coverage,
                                        extra_params=uniprot_fast_params)
        if uniprot_hit:
            uid, _, _ = uniprot_hit
            uniprot_acc = uid.split('|')[1] if '|' in uid else uid
            tmp_pdb = os.path.join(temp_dir, f"AF-{uniprot_acc}.pdb")
            if get_alphafold_structure(uniprot_acc, tmp_pdb):
                af_seq = extract_seq_from_pdb(tmp_pdb, 'A')
                if af_seq:
                    ok, *_ , start, end = align_and_verify(
                        sequence, af_seq,
                        args.min_identity, args.min_coverage)
                    if ok:
                        truncate_and_save_pdb(tmp_pdb, output_pdb,
                                             'A', start, end)
                        print("\n🎉 Done!")
                        shutil.rmtree(temp_dir)
                        sys.exit(0)

        # Step 3 – Boltz-2 de-novo prediction (DOCKER)
        boltz_dir = os.path.join(output_dir, f"{base_name}_boltz_results")
        if run_boltz(input_path, output_pdb, boltz_dir,
                     args.docker_image, args.boltz_cache):
            print("\n🎉 Done!")
            sys.exit(0)
        else:
            print("\n❌ Pipeline failed.", file=sys.stderr)
            sys.exit(1)

    try:
        shutil.rmtree(temp_dir)
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
