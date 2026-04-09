#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structure Predictor – Docker wrapper (v2)

Thin host-side script that runs the full prediction pipeline inside a
Docker container with GPU support.  No Conda environment or Python
scientific packages are needed on the host – only Python 3 and Docker.

Usage:
    python3 structure_predictor_docker.py <input_file> [options]

Requirements on the host:
    • Docker >= 19.03
    • NVIDIA driver >= 525.60.13
    • nvidia-container-toolkit (for --gpus support)
"""

import argparse
import os
import subprocess
import sys

DOCKER_IMAGE = "biochorl/structure-predictor:latest"


def resolve(path):
    """Return the absolute, realpath-resolved version of *path*."""
    return os.path.realpath(os.path.abspath(path))


def build_docker_command(args):
    """Translate host paths into container mount-points and return the
    full ``docker run`` command as a list of strings."""

    # ── Resolve host paths ───────────────────────────────────────────────
    input_path = resolve(args.input_file)
    if not os.path.isfile(input_path):
        sys.exit(f"ERROR: Input file not found: {input_path}")

    input_dir   = os.path.dirname(input_path)
    input_fname = os.path.basename(input_path)

    pdb_db_path     = resolve(args.pdb_db)
    uniprot_db_path = resolve(args.uniprot_db)

    pdb_db_dir     = os.path.dirname(pdb_db_path)
    pdb_db_name    = os.path.basename(pdb_db_path)
    uniprot_db_dir = os.path.dirname(uniprot_db_path)
    uniprot_db_name = os.path.basename(uniprot_db_path)

    boltz_cache = resolve(args.boltz_cache)
    os.makedirs(boltz_cache, exist_ok=True)

    # ── Container mount-points ───────────────────────────────────────────
    #   /workspace       ← input/output directory  (rw)
    #   /db/pdb          ← PDB database directory  (ro)
    #   /db/uniprot      ← UniProt database dir     (ro)
    #   /root/.cache     ← Boltz model weights      (rw, persistent)
    volumes = [
        f"{input_dir}:/workspace",
        f"{pdb_db_dir}:/db/pdb:ro",
        f"{uniprot_db_dir}:/db/uniprot:ro",
        f"{boltz_cache}:/root/.cache",
    ]

    # ── Assemble inner command (arguments to structure_predictor.py) ─────
    inner_args = [
        f"/workspace/{input_fname}",
        "--min_identity", str(args.min_identity),
        "--min_coverage", str(args.min_coverage),
        "--pdb_db",      f"/db/pdb/{pdb_db_name}",
        "--uniprot_db",  f"/db/uniprot/{uniprot_db_name}",
    ]

    # ── Build docker run command ─────────────────────────────────────────
    cmd = ["docker", "run", "--rm"]

    # GPU access
    cmd += ["--gpus", "all"]

    # Volume mounts
    for v in volumes:
        cmd += ["-v", v]

    # Working directory inside container
    cmd += ["-w", "/workspace"]

    # Image
    cmd.append(args.docker_image)

    # Script arguments
    cmd += inner_args

    return cmd


def main():
    parser = argparse.ArgumentParser(
        description="Structure Predictor – Docker wrapper (v2). "
                    "Runs the prediction pipeline inside a GPU-capable "
                    "Docker container.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic run
  python3 structure_predictor_docker.py my_protein.fasta \\
      --pdb_db /data/mmseqs_dbs/PDB \\
      --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL

  # Custom thresholds
  python3 structure_predictor_docker.py my_protein.fasta \\
      --min_identity 95.0 --min_coverage 90.0 \\
      --pdb_db /data/mmseqs_dbs/PDB \\
      --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL

  # Build the Docker image first
  python3 structure_predictor_docker.py my_protein.fasta --build \\
      --pdb_db /data/mmseqs_dbs/PDB \\
      --uniprot_db /data/mmseqs_dbs/UniProtKB-TrEMBL
""")

    parser.add_argument("input_file",
                        help="Input file (.fasta, .fa, .pdb, or .cif).")
    parser.add_argument("--min_identity", type=float, default=99.0,
                        help="Minimum sequence identity (%%) [default: 99.0].")
    parser.add_argument("--min_coverage", type=float, default=99.0,
                        help="Minimum query coverage (%%) [default: 99.0].")
    parser.add_argument("--pdb_db", type=str, required=True,
                        help="Absolute path to the local MMseqs2-formatted "
                             "PDB database (e.g. /data/mmseqs_dbs/PDB).")
    parser.add_argument("--uniprot_db", type=str, required=True,
                        help="Absolute path to the local MMseqs2-formatted "
                             "UniProtKB-TrEMBL database.")
    parser.add_argument("--docker_image", type=str, default=DOCKER_IMAGE,
                        help=f"Docker image to use [default: {DOCKER_IMAGE}].")
    parser.add_argument("--boltz_cache", type=str,
                        default=os.path.expanduser("~/.cache/boltz_docker"),
                        help="Host directory for persistent Boltz model "
                             "weight cache [default: ~/.cache/boltz_docker].")
    parser.add_argument("--build", action="store_true",
                        help="Build (or rebuild) the Docker image before "
                             "running the pipeline.")

    args = parser.parse_args()

    # ── Optional build step ──────────────────────────────────────────────
    if args.build:
        script_dir = os.path.dirname(resolve(__file__))
        dockerfile = os.path.join(script_dir, "Dockerfile")
        if not os.path.isfile(dockerfile):
            sys.exit(f"ERROR: Dockerfile not found at {dockerfile}")
        print(f"🐳 Building Docker image '{args.docker_image}' …")
        build_cmd = [
            "docker", "build",
            "-t", args.docker_image,
            "-f", dockerfile,
            script_dir,
        ]
        result = subprocess.run(build_cmd)
        if result.returncode != 0:
            sys.exit("ERROR: Docker build failed.")
        print(f"✔️  Image '{args.docker_image}' built successfully.\n")

    # ── Run ──────────────────────────────────────────────────────────────
    cmd = build_docker_command(args)

    print("🐳 Launching Docker container …")
    print(f"   Image : {args.docker_image}")
    print(f"   Input : {resolve(args.input_file)}")
    print(f"   PDB DB: {resolve(args.pdb_db)}")
    print(f"   UniProt DB: {resolve(args.uniprot_db)}")
    print(f"   Cache : {resolve(args.boltz_cache)}")
    print()

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
