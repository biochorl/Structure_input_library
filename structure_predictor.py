#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import shutil
import requests
import json
import time
import datetime
import gzip
import re
from glob import glob

try:
    from Bio import SeqIO, Align
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.MMCIFParser import MMCIFParser
except ImportError:
    print("ERRORE: Biopython non è installato. Per favore, installalo con 'pip install biopython'", file=sys.stderr)
    sys.exit(1)

CONDA_ENV_PATH = "/home/marco/miniconda3/envs/Boltz-2"
RESIDUE_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", 
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", 
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", 
    "TYR": "Y", "VAL": "V"
}

class Logger:
    def __init__(self, log_file_path):
        self.terminal = sys.stdout
        self.log_file = open(log_file_path, 'w') # Open in write mode to overwrite old logs

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

def extract_seq_from_pdb(pdb_path, chain_id):
    print(f"🧬 Estraggo la sequenza dalla catena '{chain_id}' del file PDB/CIF...")
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
        seq = [RESIDUE_MAP.get(r.get_resname().strip().upper(), 'X') for r in chain.get_residues() if r.get_id()[0] == ' ']
        return "".join(seq)
    except Exception as e: 
        print(f"ERRORE estrazione sequenza PDB: {e}")
        return None

def align_and_verify(query_seq, target_seq, min_identity, min_coverage):
    print("\n" + "="*35)
    print("🔬  INIZIO ALLINEAMENTO DIRETTO  🔬")
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
    print("\n[Visualizzazione Miglior Allineamento]:\n" + str(best_alignment))

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
    print("📊 Risultati Allineamento:")
    print(f"   - Identità: {identity:.2f}%")
    print(f"   - Copertura: {coverage:.2f}%")
    print(f"   - Regione allineata nel PDB: Residui {t_start_pdb} - {t_end_pdb}")
    print("-"*35 + "\n")

    is_match = identity >= min_identity and coverage >= min_coverage
    if is_match: return True, identity, coverage, t_start_pdb, t_end_pdb
    return False, identity, coverage, None, None

def truncate_and_save_pdb(input_pdb_path, output_pdb_path, chain_id, start, end):
    print(f"✂️  Troncatura del modello (Residui: {start}-{end}, Catena: {chain_id})...")
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
        print(f"✔️ Modello troncato salvato: {output_pdb_path}")
        return True
    except Exception as e: 
        print(f"ERRORE durante la troncatura: {e}", file=sys.stderr)
        shutil.copy(input_pdb_path, output_pdb_path)
        return False

def process_pdb_file(pdb_path, output_path):
    print(f"📖 Elaborazione del file di struttura: {pdb_path}")
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
        print(f"✔️ Catena '{chain.id}' salvata in: {output_path}")
        return True
    except Exception as e: print(f"ERRORE elaborazione PDB/CIF: {e}", file=sys.stderr); return False

def check_mmseqs_availability():
    if not shutil.which("mmseqs"): sys.exit("ERRORE: 'mmseqs' non trovato. Installare MMseqs2.")
    print("✔️ Comando 'mmseqs' trovato.")

def parse_fasta(fasta_path):
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq).upper(), len(record.seq)
    except Exception as e: sys.exit(f"ERRORE lettura FASTA: {e}")

def run_mmseqs_search(temp_fasta, query_len, db_path, temp_dir, min_identity, min_coverage, extra_params=None):
    db_name = os.path.basename(db_path)
    print(f"\n🔎 Avvio ricerca MMseqs2 su '{db_name}'...")
    out_m8 = os.path.join(temp_dir, f"search_{db_name}.m8")
    tmp_mmseqs = os.path.join(temp_dir, f"mmseqs_tmp_{db_name}")
    os.makedirs(tmp_mmseqs, exist_ok=True)
    command = ["mmseqs", "easy-search", temp_fasta, db_path, out_m8, tmp_mmseqs, "--format-output", "target,pident,alnlen,qlen,tlen,tstart,tend"]
    if extra_params:
        command.extend(extra_params)
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        if not os.path.exists(out_m8): print("MMseqs2 non ha prodotto risultati."); return None
        with open(out_m8, 'r') as f: hits = f.read().strip().split('\n')
        if not hits or not hits[0]: print("MMseqs2 non ha prodotto risultati."); return None
        for line in hits:
            parts = line.split('\t')
            if len(parts) < 7: continue
            sseqid, pident, align_len, qlen, slen, sstart, send = parts
            coverage = (int(align_len) / query_len) * 100
            if float(pident) >= min_identity and coverage >= min_coverage:
                print(f"✔️ Trovata corrispondenza in {db_name}: {sseqid.strip()}")
                return (sseqid.strip(), int(sstart), int(send))
        print(f"\n❌ Nessuna corrispondenza in {db_name} ha superato le soglie.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"ERRORE MMseqs2 (exit code {e.returncode})", file=sys.stderr)
        if e.stderr: print(f"  stderr: {e.stderr.strip()}", file=sys.stderr)
        return None
    except Exception as e: print(f"ERRORE MMseqs2: {e}", file=sys.stderr); return None

def download_pdb_file(pdb_id, output_path, source="rcsb"):
    url = {"rcsb": f"https://files.rcsb.org/download/{pdb_id}.pdb", "alphafold": f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb"}.get(source)
    print(f"⬇️  Download del file PDB da {url}...")
    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(output_path, 'wb') as f: shutil.copyfileobj(r.raw, f)
            print(f"✔️ File salvato con successo in: {output_path}")
            return True
    except requests.exceptions.RequestException as e: print(f"ERRORE: Download fallito. {e}", file=sys.stderr); return False

def get_alphafold_structure(uniprot_id, output_path):
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"\n🔎 Verifica del modello pre-calcolato su AlphaFold DB (API): {api_url}...")
    try:
        res = requests.get(api_url, timeout=60)
        if res.status_code != 200:
            print(f"❌ Richiesta API fallita con codice di stato: {res.status_code}", file=sys.stderr)
            print(f"   Contenuto della risposta: {res.text}", file=sys.stderr)
            return False
        try:
            prediction_data = res.json()
        except requests.exceptions.JSONDecodeError:
            print(f"❌ La risposta dell'API non è un JSON valido.", file=sys.stderr)
            print(f"   Contenuto della risposta: {res.text}", file=sys.stderr)
            return False
        if not prediction_data:
            print("❌ L'API ha restituito una risposta vuota.", file=sys.stderr)
            return False
        entry_data = prediction_data[0]
        model_url = entry_data.get('pdbUrl') or entry_data.get('cifUrl')
        avg_plddt = entry_data.get('globalMetricValue')
        if not model_url or avg_plddt is None:
            print("❌ La risposta dell'API non contiene un URL per il modello o il punteggio plddt.", file=sys.stderr)
            print(f"   Dati ricevuti: {entry_data}", file=sys.stderr)
            return False
        print(f"📊 Controllo qualità AlphaFold DB: pLDDT={avg_plddt:.2f}")
        if avg_plddt >= 75:
            print("✔️ Il modello soddisfa i criteri di qualità (pLDDT >= 75).")
            print(f"⬇️  Download del modello da {model_url}...")
            with requests.get(model_url, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(output_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)
                print(f"✔️ File salvato con successo in: {output_path}")
                return True
        else:
            print("❌ Il modello non supera la soglia di qualità.")
            return False
    except Exception as e:
        print(f"ERRORE imprevisto durante la chiamata API ad AlphaFold DB: {e}", file=sys.stderr)
        return False

def run_boltz(fasta_path, final_output_path, boltz_output_dir):
    print(f"\n🚀 Avvio della predizione con Boltz...")
    if os.path.exists(boltz_output_dir): shutil.rmtree(boltz_output_dir)
    os.makedirs(boltz_output_dir, exist_ok=True)
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        temp_fasta_path = os.path.join(boltz_output_dir, "temp_input.fasta")
        with open(temp_fasta_path, "w") as temp_f: temp_f.write(f">A|protein|\n{record.seq}\n")
        command = ["boltz", "predict", os.path.abspath(temp_fasta_path), "--use_potentials", "--diffusion_samples", "5", "--use_msa_server"]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=3600, cwd=boltz_output_dir)
        prediction_base_dir = os.path.join(boltz_output_dir, f"boltz_results_{os.path.splitext(os.path.basename(temp_fasta_path))[0]}", "predictions", os.path.splitext(os.path.basename(temp_fasta_path))[0])
        score_files = glob(os.path.join(prediction_base_dir, "confidence_*.json"))
        ranked_models = sorted([(json.load(open(f)).get('confidence_score'), f.replace("confidence_", "").replace(".json", ".cif")) for f in score_files], reverse=True)
        best_score, best_model_file = ranked_models[0]
        print(f"\n📊 Miglior punteggio di confidenza: {best_score:.4f}")
        if best_score >= 0.8:
            print("✔️ Il modello soddisfa i criteri di qualità.")
            return process_pdb_file(best_model_file, final_output_path)
        else: print("❌ QUALITÀ INSUFFICIENTE."); return False
    except Exception as e: print(f"ERRORE Boltz: {e}", file=sys.stderr); return False

def main():
    conda_env_path = os.path.abspath(CONDA_ENV_PATH)
    current_conda_env = os.environ.get("CONDA_PREFIX")
    if current_conda_env: current_conda_env = os.path.abspath(current_conda_env)

    if current_conda_env != conda_env_path:
        print(f"🐍 Ambiente Conda '{os.path.basename(conda_env_path)}' non attivo.", flush=True)
        print("Tentativo di riavvio automatico dello script nell'ambiente corretto...", flush=True)
        command = ["conda", "run", "-p", conda_env_path, sys.executable, *sys.argv]
        try:
            result = subprocess.run(command, check=True)
            sys.exit(result.returncode)
        except (FileNotFoundError, subprocess.CalledProcessError, Exception) as e:
            print(f"ERRORE critico durante il tentativo di attivazione dell'ambiente: {e}", file=sys.stderr)
            sys.exit(1)

    parser = argparse.ArgumentParser(description="Processa o predice strutture proteiche.")
    parser.add_argument("input_file", help="File di input (.fasta, .pdb, .cif).")
    parser.add_argument("--min_identity", type=float, default=99.0, help="Identità minima (%%).")
    parser.add_argument("--min_coverage", type=float, default=99.0, help="Copertura minima (%%).")
    parser.add_argument("--pdb_db", type=str, default="/media/marco/Database/Test_Database/PDB", help="Percorso al database PDB per MMseqs2.")
    parser.add_argument("--uniprot_db", type=str, default="/media/marco/Database/Test_Database/UniProtKB-TrEMBL", help="Percorso al database UniProtKB per MMseqs2.")
    args = parser.parse_args()

    input_path = os.path.abspath(args.input_file)
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    output_dir = os.path.dirname(input_path)
    output_pdb_path = os.path.join(output_dir, f"{base_name}.pdb")
    log_file = os.path.join(output_dir, f"{base_name}_report.log")
    
    sys.stdout = Logger(log_file)
    
    print(f"✔️ Ambiente Conda corretto ('{os.path.basename(conda_env_path)}') attivo.")

    temp_dir = os.path.join(output_dir, f"{base_name}_temp_files")
    if os.path.exists(temp_dir): shutil.rmtree(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)

    print(f"Avvio elaborazione per: {os.path.basename(input_path)}")

    if input_path.endswith(('.pdb', '.cif')):
        process_pdb_file(input_path, os.path.join(output_dir, f"{base_name}_chain1.pdb"))
    elif input_path.endswith(('.fasta', '.fa')):
        check_mmseqs_availability()
        sequence, seq_len = parse_fasta(input_path)
        print(f"Sequenza letta: {seq_len} aa.")

        temp_query_fasta = os.path.join(temp_dir, "query.fasta")
        with open(temp_query_fasta, "w") as f_out: f_out.write(f">query\n{sequence}\n")

        pdb_db = os.path.abspath(args.pdb_db)
        uniprot_db = os.path.abspath(args.uniprot_db)

        if (pdb_hit := run_mmseqs_search(temp_query_fasta, seq_len, pdb_db, temp_dir, args.min_identity, args.min_coverage)):
            sseqid, _, _ = pdb_hit
            parts = sseqid.replace('_', '|').split('|')
            if len(parts) >= 2:
                # Gestisce sia pdb|1ABC|A che 1ABC_A
                pdb_code, chain_id = parts[-2][:4].lower(), parts[-1]
                if len(pdb_code) == 4:
                    temp_pdb_path = os.path.join(temp_dir, f"{pdb_code}.pdb")
                    if download_pdb_file(pdb_code, temp_pdb_path, source="rcsb"):
                        pdb_seq = extract_seq_from_pdb(temp_pdb_path, chain_id)
                        if pdb_seq:
                            is_match, identity, coverage, start, end = align_and_verify(sequence, pdb_seq, args.min_identity, args.min_coverage)
                            if is_match:
                                truncate_and_save_pdb(temp_pdb_path, output_pdb_path, chain_id, start, end)
                                print("\n🎉 Operazione completata!"); shutil.rmtree(temp_dir); sys.exit(0)
        
        # Use faster search parameters for the large UniProt database:
        # -s 4 reduces sensitivity (default 5.7) for ~3x speedup,
        # --max-seqs 100 limits prefilter results since we only need the top hit.
        uniprot_fast_params = ["-s", "4", "--max-seqs", "100"]
        if (uniprot_hit := run_mmseqs_search(temp_query_fasta, seq_len, uniprot_db, temp_dir, args.min_identity, args.min_coverage, extra_params=uniprot_fast_params)):
            uniprot_id, _, _ = uniprot_hit
            # Gestisce formati come sp|P12345|... o tr|A0A0|... o semplici ID
            uniprot_acc = uniprot_id.split('|')[1] if '|' in uniprot_id else uniprot_id
            temp_pdb_path = os.path.join(temp_dir, f"AF-{uniprot_acc}.pdb")
            if get_alphafold_structure(uniprot_acc, temp_pdb_path):
                alphafold_seq = extract_seq_from_pdb(temp_pdb_path, 'A')
                if alphafold_seq:
                    is_match, identity, coverage, start, end = align_and_verify(sequence, alphafold_seq, args.min_identity, args.min_coverage)
                    if is_match:
                        truncate_and_save_pdb(temp_pdb_path, output_pdb_path, 'A', start, end)
                        print("\n🎉 Operazione completata!"); shutil.rmtree(temp_dir); sys.exit(0)


        boltz_results_dir = os.path.join(output_dir, f"{base_name}_boltz_results")
        if run_boltz(input_path, output_pdb_path, boltz_results_dir):
            print("\n🎉 Operazione completata!")
            sys.exit(0)
        else:
            print("\n❌ Operazione fallita.", file=sys.stderr)
            sys.exit(1)
    
    try:
        shutil.rmtree(temp_dir)
    except FileNotFoundError:
        pass

if __name__ == "__main__":
    main()
