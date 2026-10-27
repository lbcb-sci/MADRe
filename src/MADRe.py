#!/usr/bin/env python3

import argparse
import configparser
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import threading
import time
from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor
from NovelPlacement import ancestor, read_taxonomy

# MADRe pipeline:
#   1 assembly (metaFlye / metaMDBG / Myloasm), or given contigs (--assembly)
#   2 assembly contigs mapped to the full reference database
#   3 collapsed strains per contig (HairSplitter)
#   4 contig classification (KNOWN / UNCERTAIN / NOVEL_CANDIDATE) and the
#     reduced database, selected as whole genomes (DatabaseReduction,
#     needs --reference_genome_map)
#   5 grouping of UNCERTAIN / novel contigs (NovelDetection)
#   6 extended database = reduced database + UNCERTAIN / novel contigs
#     (+ host reference)
#   7 reads mapped to the extended database
#   8 StrainRefine: novel-contig placement in the taxonomy (on by
#     default, needs --nodes/--names; --no-novel-placement turns it off), read
#     assignment, strain clustering
#   9 CalculateAbundances: abundance, species summary and taxonomic profile
# --min-abundance leaves species below it out of the final taxonomic profile.
# --nodes/--names are used by steps 8 and 9 only; step 5 runs without taxonomy.
import CalculateAbundances
import DatabaseReduction as DR
import NovelDetection
import StrainRefine


# =============================================================================
# Defaults
# =============================================================================

# Shown when --reference_genome_map is missing.
REFERENCE_GENOME_MAP_EXAMPLE = """\
Expected: a tab-separated file with a header line and one row per sequence of
the reference database. genome_id groups the sequences of one genome (e.g. an
assembly accession); genome_length is that whole genome's length in bp (sum
of all its sequences), repeated on each of its rows. Example:

reference\tgenome_id\tgenome_length
|1000373|NC_016757.1\tGCF_000895895.1\t17078
|1000373|NC_016758.1\tGCF_000895895.1\t17078
|562|NZ_CM137271.1\tGCF_049197725.1\t5411416

The reference column must match the database FASTA headers exactly
(the text after '>' up to the first space)."""

MIN_CONTIG_LEN = 1000
MAX_COLLAPSED_STRAINS_OVERHEAD = 2

DEFAULT_KNOWN_THR = 0.75
DEFAULT_NOVEL_THR = 0.35

DEFAULT_CLUSTER_EPS = 0.8
DEFAULT_HIGH_SCORE_THRESHOLD = 0.6
DEFAULT_MIN_HIGH_SCORE_READS = 5
DEFAULT_MID_HIGH_SCORE_READS = 10
DEFAULT_MIN_HIGH_SCORE_FRACTION = 0.8
DEFAULT_SPECIES_MIN_READ_COUNT = 5
DEFAULT_SPECIES_MIN_MEAN_SCORE = 0.6
DEFAULT_SPECIES_LOW_COUNT_CAP = 30

# minimap2 (2.28: index.c, options.c) indexes a large database in parts: it reads the targets in mini-batches of at
# least 50 Mbp and closes a part once it holds more than -I bases (default 8G). Without --split-prefix each part is
# mapped on its own and its alignments are written before the next part is read, so mapping the parts separately
# and joining the outputs in part order gives the same PAF as one run. One run computes the repetitive-minimizer
# cutoff (mid_occ) from the first part only and keeps it for all parts, so the other parts get it with -f.
MINIMAP_INDEX_PART_BASES = 8_000_000_000
MINIMAP_INDEX_MINI_BATCH = 50_000_000
DEFAULT_DB_MAP_JOBS = 8


# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[logging.StreamHandler()])


# =============================================================================
# Helpers
# =============================================================================


def run_command(cmd, target_file=None, force=False, stderr_path=None):
    """
    Run an external command unless target_file already exists.

    cmd must be a list.
    """

    if target_file is not None and os.path.exists(target_file) and not force:
        logging.info("Skipping step because output already exists: %s", target_file)
        return

    logging.info("Running: %s", " ".join(shlex.quote(str(x)) for x in cmd))

    stderr_handle = None

    try:
        if stderr_path:
            stderr_handle = open(stderr_path, "w")

        subprocess.run(cmd, check=True, stderr=stderr_handle)

    finally:
        if stderr_handle:
            stderr_handle.close()

    if target_file is not None and not os.path.exists(target_file):
        raise RuntimeError("Command completed but expected output " f"was not created: {target_file}")


def run_command_stdout(cmd, output_path, force=False, stderr_path=None):
    """
    Run a command and redirect stdout to output_path.

    Output is written to output_path + ".tmp" and renamed to output_path only
    after the command finished successfully, so an interrupted or failed run
    never leaves a partial file that a later run would take as finished.
    """

    if os.path.exists(output_path) and not force:
        logging.info("Skipping step because output already exists: %s", output_path)
        return

    logging.info("Running: %s > %s", " ".join(shlex.quote(str(x)) for x in cmd), output_path)

    stderr_handle = None
    tmp_path = str(output_path) + ".tmp"

    try:
        if stderr_path:
            stderr_handle = open(stderr_path, "w")

        with open(tmp_path, "w") as stdout_handle:

            subprocess.run(cmd, stdout=stdout_handle, stderr=stderr_handle, check=True)

        os.replace(tmp_path, output_path)

    finally:
        if stderr_handle:
            stderr_handle.close()
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    if not os.path.exists(output_path):
        raise RuntimeError(f"Expected output was not created: {output_path}")


def database_index_parts(database):
    """
    The parts minimap2 splits database's index into, from database.fai: a list of (first sequence name, start byte,
    size in bytes or None for the rest of the file, number of sequences). None when the database cannot be read
    by byte range (no .fai, or a .gz without the .gzi of bgzip, or bgzip missing).
    """

    if not os.path.exists(database + ".fai"):
        return None
    if str(database).endswith(".gz") and (not os.path.exists(database + ".gzi") or not shutil.which("bgzip")):
        return None

    parts = []
    first_name, start, n_seq, part_bases, mini_bases = None, 0, 0, 0, 0
    with open(database + ".fai") as fai:
        for line in fai:
            name, length, offset, line_bases, line_width = line.split("\t")[:5]
            length, offset, line_bases, line_width = int(length), int(offset), int(line_bases), int(line_width)
            # end of the record: full lines with their newline, then the last partial line with its newline
            full_lines, rest = divmod(length, line_bases) if line_bases else (0, 0)
            end = offset + full_lines * line_width + (rest + line_width - line_bases if rest else 0)
            if first_name is None:
                first_name = name
            n_seq += 1
            part_bases += length
            mini_bases += length
            if mini_bases >= MINIMAP_INDEX_MINI_BATCH:
                mini_bases = 0
                if part_bases > MINIMAP_INDEX_PART_BASES:
                    parts.append((first_name, start, end - start, n_seq))
                    first_name, start, n_seq, part_bases = None, end, 0, 0
    if n_seq:
        parts.append((first_name, start, None, n_seq))
    return parts


def database_part_reader(database, start, size):
    """Shell command writing size bytes of database from byte start (to the end when size is None) to stdout."""

    db = shlex.quote(str(database))
    if str(database).endswith(".gz"):
        return f"bgzip -c -d -b {start}" + (f" -s {size}" if size is not None else "") + f" {db}"
    return f"dd if={db} bs=4M iflag=skip_bytes,count_bytes skip={start}" + (f" count={size}" if size is not None else "") + " status=none"


def map_assembly_to_database(minimap, database, assembly, output_path, threads, jobs, force=False, stderr_path=None):
    """
    Step 2: minimap2 -x asm5 of the assembly against the full database.

    One minimap2 run spends most of its time building the index parts one after another (reading and sketching the
    database on a few cores; mapping the contigs to a part takes seconds). Here the parts (database_index_parts) are
    streamed from the database by byte range and mapped by up to `jobs` minimap2 runs at a time, and their PAFs are
    written in part order: the same alignments, in the same order, as one run. The first part starts alone; once its
    index is built, its mid_occ (from its log) is given to the other parts with -f, as one run would use it.
    No files besides the PAF (written to <output>.tmp and renamed, as in run_command_stdout) and the log: a part's
    alignments and log are held in memory until the parts before it are written. Falls back to one run when the
    database has no byte-range access (see database_index_parts) or only one part.
    """

    if os.path.exists(output_path) and not force:
        logging.info("Skipping step because output already exists: %s", output_path)
        return

    parts = database_index_parts(database) if jobs > 1 else None
    if not parts or len(parts) == 1:
        if jobs > 1 and parts is None:
            logging.info("Database cannot be read by byte range (needs .fai, and .gzi + bgzip for .gz); one minimap2 run.")
        cmd = [minimap, "-x", "asm5", database, assembly, "-t", str(threads)]
        run_command_stdout(cmd, output_path, force=force, stderr_path=stderr_path)
        return

    jobs = min(jobs, len(parts))
    part_threads = max(1, threads // jobs)
    logging.info("Mapping to %d database index parts, %d at a time with %d threads each.", len(parts), jobs, part_threads)

    # every part must start with the header of its first sequence, else the .fai does not describe the file
    for i, (name, start, size, n_seq) in enumerate(parts):
        head = subprocess.run(["bash", "-o", "pipefail", "-c", database_part_reader(database, start, len(name) + 2)], capture_output=True, check=True)
        if head.stdout.split()[0] != b">" + name.encode():
            raise RuntimeError(f"Database part {i + 1} does not start with >{name} at byte {start}; is {database}.fai up to date?")

    first_mid_occ = []
    first_mid_occ_found = threading.Event()

    def map_part(i, mid_occ=None):
        """(PAF bytes, log bytes) of part i."""
        name, start, size, n_seq = parts[i]
        began = time.time()
        cmd = (
            f"{database_part_reader(database, start, size)} | "
            f"{shlex.quote(minimap)} -x asm5 -I {MINIMAP_INDEX_PART_BASES}"
            + (f" -f {mid_occ}" if mid_occ is not None else "")
            + f" -t {part_threads} - {shlex.quote(str(assembly))}"
        )
        proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_lines = []

        def read_log():
            for line in proc.stderr:
                log_lines.append(line)
                found = re.search(rb"mid_occ = (\d+)", line) if i == 0 else None
                if found:
                    first_mid_occ.append(int(found.group(1)))
                    first_mid_occ_found.set()

        log_reader = threading.Thread(target=read_log)
        log_reader.start()
        paf = proc.stdout.read()
        proc.wait()
        log_reader.join()
        if i == 0:
            first_mid_occ_found.set()  # also when it failed before printing mid_occ
        if proc.returncode != 0:
            raise RuntimeError(f"minimap2 failed on database part {i + 1} (exit {proc.returncode}): {cmd}\n" + b"".join(log_lines[-20:]).decode(errors="replace"))
        logging.info("Database part %d/%d mapped (%d sequences from %s, %.0f s).", i + 1, len(parts), n_seq, name, time.time() - began)
        return paf, b"".join(log_lines)

    tmp_path = str(output_path) + ".tmp"
    log_handle = open(stderr_path, "wb") if stderr_path else None
    try:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            results = [pool.submit(map_part, 0)]
            first_mid_occ_found.wait()
            if not first_mid_occ:
                results[0].result()  # raises if the first part failed
                raise RuntimeError("No mid_occ in the log of database part 1; cannot map the other parts the way one minimap2 run would.")
            logging.info("mid_occ of the first database part: %d (used for all parts).", first_mid_occ[0])
            results += [pool.submit(map_part, i, first_mid_occ[0]) for i in range(1, len(parts))]
            with open(tmp_path, "wb") as out:
                for i, result in enumerate(results):
                    try:
                        paf, log = result.result()
                    except Exception:
                        for waiting in results[i + 1 :]:
                            waiting.cancel()  # parts not started yet; running ones finish first
                        raise
                    out.write(paf)
                    if log_handle:
                        name, start, size, n_seq = parts[i]
                        log_handle.write(f"### database index part {i + 1}/{len(parts)}: {n_seq} sequences from {name}\n".encode())
                        log_handle.write(log)
                    results[i] = None  # free the part's alignments
        os.replace(tmp_path, output_path)
    finally:
        if log_handle:
            log_handle.close()
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def concatenate_fastas(input_fastas, output_path, force=False):
    """
    Concatenate uncompressed FASTA files (written to output_path + ".tmp"
    and renamed when complete, so an interrupted run leaves no partial file).
    """

    if os.path.exists(output_path) and not force:
        logging.info("Skipping extended DB construction because " "output already exists: %s", output_path)
        return

    logging.info("Building extended database: %s", output_path)

    tmp_path = str(output_path) + ".tmp"

    with open(tmp_path, "w") as out:

        for fasta in input_fastas:

            require_file(fasta, "Input FASTA")

            logging.info("  Adding %s", fasta)

            with open(fasta) as inp:
                for line in inp:
                    out.write(line)

    os.replace(tmp_path, output_path)

    if not os.path.exists(output_path):
        raise RuntimeError(f"Failed to create extended database: {output_path}")


HOST_PREFIX = "HOST_"


def prepare_host_reference(host_fasta, host_name, out_folder, force=False):
    """
    Write the host reference with MADRe host headers:

        >HOST|HOST_<name>|<original id> <original description>

    StrainRefine uses the field after the first '|' as the reference id, so
    every host sequence shares the id HOST_<name>, which is registered as its
    own group in the extended strain-species JSON (reads assigned to it are
    counted as host, never as a taxon). Accepts plain or gzipped FASTA.
    Returns (path to rewritten FASTA, host group id).
    """
    require_file(host_fasta, "Host reference FASTA")
    safe = "".join(c if c.isalnum() else "_" for c in host_name).strip("_") or "host"
    group_id = HOST_PREFIX + safe
    out_path = os.path.join(out_folder, f"host_reference.{safe}.fa")
    if os.path.exists(out_path) and not force:
        logging.info("Host reference already prepared: %s", out_path)
        return out_path, group_id

    import gzip

    opener = gzip.open if str(host_fasta).endswith(".gz") else open
    n_seq = 0
    tmp_path = out_path + ".tmp"
    with opener(host_fasta, "rt") as inp, open(tmp_path, "w") as out:
        for line in inp:
            if line.startswith(">"):
                seq_id, _, desc = line[1:].rstrip("\n").partition(" ")
                seq_id = seq_id.replace("|", "_")
                out.write(f">HOST|{group_id}|{seq_id}" + (f" {desc}" if desc else "") + "\n")
                n_seq += 1
            else:
                out.write(line)
    if n_seq == 0:
        os.remove(tmp_path)
        raise RuntimeError(f"No FASTA records found in host reference: {host_fasta}")
    os.replace(tmp_path, out_path)
    logging.info("Host reference: %d sequences from %s -> %s (group %s)", n_seq, host_fasta, out_path, group_id)
    return out_path, group_id


def file_sha256(path):
    """SHA-256 of a file's content (to tell whether an input changed between runs)."""
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path, description):
    """
    Check that a required file exists.
    """

    if path is None:
        raise FileNotFoundError(f"{description} was not provided.")

    if not os.path.exists(path):
        raise FileNotFoundError(f"{description} does not exist: {path}")


# =============================================================================
# Assembly
# =============================================================================


METAMDBG_GRAPH_K = 21


def metamdbg_graph(path_metamdbg, assembly_dir, threads, force, log_path):
    """
    Build metaMDBG's assembly graph ("metaMDBG gfa") for contig links in step 5.

    Uses k = METAMDBG_GRAPH_K (every contig has a path there), or the smallest available k.
    Keeps the graph without sequences and the contig -> unitig table; the graph with
    sequences is removed. Returns (gfa, contig_paths), or (None, None) if it fails.
    """
    listing = subprocess.run([path_metamdbg, "gfa", "--assembly-dir", assembly_dir, "--k", "0"], capture_output=True, text=True)
    available = [int(k) for k in re.findall(r"^\s*-\s*(\d+)\s*\(", listing.stdout + listing.stderr, flags=re.M)]
    if not available:
        logging.warning("metaMDBG gfa: no k values available; step 5 runs without assembly-graph links.")
        return None, None
    k = METAMDBG_GRAPH_K if METAMDBG_GRAPH_K in available else min(available)
    gfa = os.path.join(assembly_dir, f"assemblyGraph_k{k}.noseq.gfa")
    contig_paths = os.path.join(assembly_dir, f"assemblyGraph_k{k}_contigPath.tsv")
    if force or not (os.path.exists(gfa) and os.path.exists(contig_paths)):
        logging.info("Building metaMDBG assembly graph (k=%d)", k)
        with open(log_path, "a") as log:
            result = subprocess.run([path_metamdbg, "gfa", "--assembly-dir", assembly_dir, "--k", str(k), "--threads", str(threads)], stdout=log, stderr=log)
        if result.returncode != 0 or not (os.path.exists(gfa) and os.path.exists(contig_paths)):
            logging.warning("metaMDBG gfa failed; step 5 runs without assembly-graph links (see %s).", log_path)
            return None, None
        with_sequences = os.path.join(assembly_dir, f"assemblyGraph_k{k}.gfa")
        if os.path.exists(with_sequences):
            os.remove(with_sequences)
    return gfa, contig_paths


def run_assembly(args, path_metaflye, path_metamdbg, path_myloasm):
    """
    Run assembly.

    With --assembly, the given contigs are used instead (copied, or decompressed, to
    <out>/external_assembly/assembly.fasta); there is no assembly graph or coverage, so
    novel-contig grouping (step 5) runs without graph links.

    Returns
    -------
    dict
        output_dir
        assembly
        assembler      metaflye / myloasm / metamdbg (None for --assembly)
        assembly_info  metaFlye assembly_info.txt
        gfa            assembly graph (metaFlye, Myloasm; metaMDBG without sequences)
        contig_paths   metaMDBG contig -> unitig table
    """

    reads_tech_flye = {"pacbio": "--pacbio-raw", "hifi": "--pacbio-hifi", "ont": "--nano-raw"}

    assembly_info = None
    assembly_gfa = None
    contig_paths = None

    if args.assembly:
        assembly_out_dir = os.path.join(args.out_folder, "external_assembly")
        assembly_fasta = os.path.join(assembly_out_dir, "assembly.fasta")
        if not os.path.exists(assembly_fasta) or args.force:
            import gzip
            import shutil

            require_file(args.assembly, "Given assembly (--assembly)")
            os.makedirs(assembly_out_dir, exist_ok=True)
            opener = gzip.open if args.assembly.endswith(".gz") else open
            with opener(args.assembly, "rb") as fin, open(assembly_fasta + ".tmp", "wb") as fout:
                shutil.copyfileobj(fin, fout)
            os.replace(assembly_fasta + ".tmp", assembly_fasta)
            logging.info("Using the given contigs instead of an assembly: %s -> %s", args.assembly, assembly_fasta)
        return {"output_dir": assembly_out_dir, "assembly": assembly_fasta, "assembler": None, "assembly_info": None, "gfa": None, "contig_paths": None}

    # Output folder, log name and assembler input option per read type.
    myloasm_setup = {"ont": ("myloasm_ont", []), "hifi": ("myloasm_hifi", ["--hifi"])}
    metamdbg_setup = {"ont": ("metaMDBG_ont", "--in-ont", "metaMDBG_ont.log"), "hifi": ("metaMDBG", "--in-hifi", "metaMDBG.log")}

    # =========================================================================
    # Myloasm (ONT or HiFi, only with --use-myloasm)
    # =========================================================================

    if args.use_myloasm:

        if not path_myloasm:
            raise ValueError("--use-myloasm was provided but myloasm " "is not defined in config.ini")

        if args.reads_flag not in myloasm_setup:
            raise ValueError("Myloasm currently supports ONT and HiFi " "in MADRe.")

        folder, extra = myloasm_setup[args.reads_flag]

        assembler = "myloasm"

        assembly_out_dir = os.path.join(args.out_folder, folder)

        assembly_fasta = os.path.join(assembly_out_dir, "assembly.fasta")

        if not os.path.exists(assembly_fasta) or args.force:

            os.makedirs(assembly_out_dir, exist_ok=True)

            command = (
                f"{shlex.quote(path_myloasm)} "
                f"--output-dir "
                f"{shlex.quote(assembly_out_dir)} "
                f"--threads {args.threads} " + "".join(f"{x} " for x in extra) + f"{shlex.quote(args.reads)} "
                f"2> "
                f"{shlex.quote(os.path.join(args.out_folder, 'myloasm.log'))} "
                f"&& cat "
                f"{shlex.quote(os.path.join(assembly_out_dir, 'assembly_primary.fa'))} "
                f"{shlex.quote(os.path.join(assembly_out_dir, 'alternate_assemblies', 'assembly_alternate.fa'))} "
                f"> {shlex.quote(assembly_fasta)}"
            )

            logging.info("Running Myloasm (primary + alternate contigs)")

            subprocess.run(command, shell=True, check=True)

        candidate_gfa = os.path.join(assembly_out_dir, "final_contig_graph.gfa")

        if os.path.exists(candidate_gfa):
            assembly_gfa = candidate_gfa

    # =========================================================================
    # metaMDBG: default for HiFi; for ONT only with --use-metamdbg
    # =========================================================================

    elif args.reads_flag == "hifi" or (args.use_metamdbg and args.reads_flag == "ont"):

        folder, input_option, log_name = metamdbg_setup[args.reads_flag]

        assembler = "metamdbg"

        assembly_out_dir = os.path.join(args.out_folder, folder)

        assembly_fasta = os.path.join(assembly_out_dir, "assembly.fasta")

        if not os.path.exists(assembly_fasta) or args.force:

            os.makedirs(assembly_out_dir, exist_ok=True)

            cmd = [path_metamdbg, "asm", "--out-dir", assembly_out_dir, input_option, args.reads, "--threads", str(args.threads)]

            run_command(cmd, force=True, stderr_path=os.path.join(args.out_folder, log_name))

            compressed = os.path.join(assembly_out_dir, "contigs.fasta.gz")

            require_file(compressed, "metaMDBG contigs")

            cmd = ["gzip", "-dc", compressed]

            run_command_stdout(cmd, assembly_fasta, force=True)

        assembly_gfa, contig_paths = metamdbg_graph(path_metamdbg, assembly_out_dir, args.threads, args.force, os.path.join(args.out_folder, log_name))

    # =========================================================================
    # metaFlye: default for ONT
    # =========================================================================

    elif args.reads_flag == "ont":

        assembler = "metaflye"

        assembly_out_dir = os.path.join(args.out_folder, "metaflye")

        assembly_fasta = os.path.join(assembly_out_dir, "assembly.fasta")

        os.makedirs(assembly_out_dir, exist_ok=True)

        cmd = [path_metaflye, reads_tech_flye[args.reads_flag], args.reads, "--out-dir", assembly_out_dir, "--threads", str(args.threads), "--meta"]

        run_command(cmd, target_file=assembly_fasta, force=args.force, stderr_path=os.path.join(args.out_folder, "metaflye.log"))

        candidate_info = os.path.join(assembly_out_dir, "assembly_info.txt")

        candidate_gfa = os.path.join(assembly_out_dir, "assembly_graph.gfa")

        if os.path.exists(candidate_info):
            assembly_info = candidate_info

        if os.path.exists(candidate_gfa):
            assembly_gfa = candidate_gfa

    else:

        raise ValueError("PacBio raw-read assembly is not currently " "implemented in MADRe.")

    require_file(assembly_fasta, "Assembly FASTA")

    return {
        "output_dir": assembly_out_dir,
        "assembly": assembly_fasta,
        "assembler": assembler,
        "assembly_info": assembly_info,
        "gfa": assembly_gfa,
        "contig_paths": contig_paths,
    }


# =============================================================================
# Main pipeline
# =============================================================================


def read_merged(path):
    """merged.dmp (old taxid -> current taxid); {} if the file does not exist."""
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        # lines without both taxids (e.g. an empty merged.dmp of a converted taxonomy) are skipped
        return {f[0]: f[1] for f in (line.replace("\t|", "").split() for line in handle) if len(f) >= 2}


def write_strain_species_json(reference_genome_map, nodes, path):
    """
    Species of every database taxid, from the taxonomy: {taxid: species taxid}, written to path (JSON) for steps 4-8.
    Taxids are the second '|' field of the genome map's references ('|562|NC_000913.3', 'kraken:taxid|562|...').
    A taxid at or above species rank maps to itself. A taxid missing from nodes.dmp is looked up under its new taxid in
    merged.dmp (next to nodes.dmp, if present); one found in neither maps to itself (with a warning).
    """
    taxids = set()
    with open(reference_genome_map) as handle:
        next(handle)
        for line in handle:
            parts = line.split("\t", 1)[0].split("|")
            if len(parts) >= 3 and parts[1]:
                taxids.add(parts[1])
    parents, ranks = read_taxonomy(nodes)
    merged = read_merged(os.path.join(os.path.dirname(os.path.abspath(nodes)), "merged.dmp"))
    current = {t: (t if t in parents else merged.get(t, t)) for t in taxids}
    n_merged = sum(1 for t in taxids if t not in parents and current[t] in parents)
    if n_merged:
        logging.info("%d database taxids were merged into newer taxids (merged.dmp); their species are taken from the newer taxid", n_merged)
    missing = sorted(t for t in taxids if current[t] not in parents)
    if missing:
        logging.warning("%d database taxids are not in %s and map to themselves, e.g. %s", len(missing), nodes, ", ".join(missing[:5]))
    mapping = {t: ancestor(current[t], parents, ranks, "species") or t for t in sorted(taxids)}
    tmp = f"{path}.tmp"
    with open(tmp, "w") as handle:
        json.dump(mapping, handle)
    os.replace(tmp, path)
    logging.info("Species of %d database taxids (%d species) from %s -> %s", len(mapping), len(set(mapping.values())), nodes, path)
    return str(path)


def main():

    parser = argparse.ArgumentParser(description="MADRe")

    # =========================================================================
    # General
    # =========================================================================

    parser.add_argument("--version", action="version", version="MADRe 2.0.0")

    parser.add_argument("--out-folder", required=True, help="Output directory.")

    parser.add_argument("--reads", required=True, help=("Input reads FASTQ/FASTA; gzip allowed " "where supported."))

    parser.add_argument(
        "--assembly",
        default=None,
        help="Contigs (FASTA, plain or .gz) to use instead of assembling the reads: step 1 is skipped. "
        "Without an assembly graph, novel-contig grouping runs without graph links.",
    )

    parser.add_argument("--reads_flag", default="ont", choices=["pacbio", "hifi", "ont"], help="Read technology. Default: ont")

    parser.add_argument("--threads", type=int, default=32, help="Number of threads. Default: 32")
    parser.add_argument(
        "--db-map-jobs",
        type=int,
        default=DEFAULT_DB_MAP_JOBS,
        help="Step 2 maps the contigs to minimap2's database index parts (~8 Gbp each) this many at a time, --threads "
        "shared between them; each holds one part's index (~50 GB RAM on a RefSeq-sized database). Same PAF as one "
        "minimap2 run. Needs the database's .fai (and .gzi for a bgzip .gz); otherwise, or with 1, one minimap2 run. "
        f"Default: {DEFAULT_DB_MAP_JOBS}",
    )

    parser.add_argument("-F", "--force", action="store_true", help=("Rerun external steps even if outputs " "already exist."))

    parser.add_argument("--config", default="config.ini", help="Path to config.ini.")

    parser.add_argument("--strictness", default="very-strict", choices=["less-strict", "strict", "very-strict"])

    parser.add_argument("--collapsed_strains_overhead", type=int, default=MAX_COLLAPSED_STRAINS_OVERHEAD)

    parser.add_argument("--min_contig_len", type=int, default=MIN_CONTIG_LEN)

    parser.add_argument("--known_thr", type=float, default=DEFAULT_KNOWN_THR, help=("vcig_norm >= this value -> KNOWN."))

    parser.add_argument("--novel_thr", type=float, default=DEFAULT_NOVEL_THR, help=("vcig_norm < this value -> NOVEL_CANDIDATE."))
    parser.add_argument(
        "--novelty-score",
        choices=["summed", "union"],
        default="union",
        help="Contig novelty score: summed = vcig of all alignments to a reference / contig length (overlapping alignments all "
        "counted); union = each contig base counted once, weighted by the best alignment covering it "
        "(union sets only the KNOWN/UNCERTAIN/NOVEL label; the best reference is chosen among the EM-selected references).",
    )

    # MAPQ=0 assembly-to-reference mappings are included by default (no MAPQ
    # filter); --mapq-filter excludes them. --no-mapq-filter states the default.
    parser.add_argument(
        "--no-mapq-filter",
        dest="no_mapq_filter",
        action="store_true",
        help=(
            "Include MAPQ=0 assembly-to-reference mappings during database " "reduction and novelty scoring (default; the option only makes " "it explicit)."
        ),
    )
    parser.add_argument(
        "--mapq-filter",
        dest="no_mapq_filter",
        action="store_false",
        help=("Exclude MAPQ=0 assembly-to-reference mappings during database " "reduction and novelty scoring."),
    )
    parser.set_defaults(no_mapq_filter=True)

    parser.add_argument("--use-myloasm", action="store_true")

    parser.add_argument(
        "--use-metamdbg",
        action="store_true",
        help=(
            "Assemble ONT reads with metaMDBG (--in-ont) instead of "
            "metaFlye. Only affects --reads_flag ont; metaMDBG already "
            "handles hifi unconditionally."
        ),
    )

    # =========================================================================
    # Taxonomy used by optional family grouping
    # =========================================================================

    taxonomy = parser.add_argument_group("Taxonomy")

    taxonomy.add_argument("--nodes", required=True, help="Taxonomy nodes.dmp (NCBI format, e.g. NCBI or GTDB in NCBI format), the version used for the database. Gives the species of every reference and is used for novel-contig placement.")

    taxonomy.add_argument("--names", required=True, help="Taxonomy names.dmp matching --nodes.")

    # =========================================================================
    # Novel grouping
    # =========================================================================

    novel = parser.add_argument_group("Novel organism grouping")

    novel.add_argument("--grouping-strict-overlap", type=float, default=0.80)

    novel.add_argument("--grouping-strict-identity", type=float, default=0.95)

    novel.add_argument("--grouping-relaxed-overlap", type=float, default=0.10)

    novel.add_argument("--grouping-relaxed-identity", type=float, default=0.70)

    # =========================================================================
    # StrainRefine
    # =========================================================================

    sr = parser.add_argument_group("StrainRefine")

    sr.add_argument("--cluster-eps", type=float, default=DEFAULT_CLUSTER_EPS)

    sr.add_argument("--high-score-threshold", type=float, default=DEFAULT_HIGH_SCORE_THRESHOLD)

    sr.add_argument("--min-high-score-reads", type=int, default=DEFAULT_MIN_HIGH_SCORE_READS)

    sr.add_argument("--mid-high-score-reads", type=int, default=DEFAULT_MID_HIGH_SCORE_READS)

    sr.add_argument("--min-high-score-fraction", type=float, default=DEFAULT_MIN_HIGH_SCORE_FRACTION)

    sr.add_argument("--species-min-read-count", type=int, default=DEFAULT_SPECIES_MIN_READ_COUNT)

    sr.add_argument("--species-min-mean-score", type=float, default=DEFAULT_SPECIES_MIN_MEAN_SCORE)

    sr.add_argument("--species-low-count-cap", type=int, default=DEFAULT_SPECIES_LOW_COUNT_CAP)

    parser.add_argument(
        "--placement-report",
        action="store_true",
        help="Draw the per-contig placement figures in novel_placement_report/ (one PNG per novel contig with a candidate family, "
        "and the numbers behind them in novel_contig_neighbours.tsv). Off by default: with thousands of novel contigs, drawing them "
        "takes hours. novel_placement.tsv is always written.",
    )
    # earlier option (the figures used to be drawn by default); accepted for old commands, no effect
    parser.add_argument("--no-placement-report", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--no-novel-placement",
        action="store_true",
        help="Do not place novel contigs in the taxonomy before StrainRefine groups reads "
        "by species (placement is on by default and needs --nodes/--names).",
    )
    parser.add_argument(
        "--keep-reduced-db",
        action="store_true",
        help="Keep reduced_db.fa (the selected database sequences). By default it is removed once extended_db.fa is built, "
        "since extended_db.fa contains the same sequences; reduced_db_refs.txt lists them.",
    )
    parser.add_argument(
        "--host-reference",
        default=None,
        help="Host reference FASTA (plain or .gz). Appended to the extended database "
        "with headers >HOST|HOST_<name>|<id>, so reads from the host are assigned "
        "to it and reported as host instead of being forced onto microbes or novel "
        "contigs. Host reads are excluded from microbial relative abundance.",
    )
    parser.add_argument("--host-name", default="host", help="Label used in host headers and outputs (HOST_<name>), e.g. human.")
    parser.add_argument(
        "--min-abundance",
        type=float,
        default=0.0,
        help="Leave species below this relative abundance (%%) out of strainrefine_taxonomic_profile.tsv; the complete profile is "
        "kept as strainrefine_taxonomic_profile.all_species.tsv. Default 0 (off).",
    )
    parser.add_argument(
        "--genome-depth",
        choices=("median", "mean"),
        default="median",
        help="Genome copies of a known genome: 'median' (default) = min(bases / genome length, length-weighted median depth "
        "of its sequences with reads), so reads piled on a small high-copy sequence (plasmid, repeat) do not inflate it; "
        "'mean' = bases / genome length (the rule of earlier MADRe versions).",
    )
    parser.add_argument(
        "--reference_genome_map",
        default=None,
        help="Required. Tab-separated table assigning every database sequence to its "
        "genome, with the whole genome's length (sum of all its sequences), e.g.: "
        "'reference  genome_id  genome_length' / "
        "'|1000373|NC_016757.1  GCF_000895895.1  17078' / "
        "'|1000373|NC_016758.1  GCF_000895895.1  17078'.",
    )
    args = parser.parse_args()
    if not args.reference_genome_map:
        parser.error(
            "--reference_genome_map is required: the database reduction uses it to "
            "select whole genomes, and StrainRefine uses it to pool the sequences "
            "of a genome and to normalise abundance by genome length.\n\n" + REFERENCE_GENOME_MAP_EXAMPLE
        )

    # =========================================================================
    # Configuration
    # =========================================================================

    logging.info("Starting MADRe...")

    logging.info("Configuration: %s", args.config)

    config = configparser.ConfigParser(inline_comment_prefixes=("#",))

    if not config.read(args.config):

        raise FileNotFoundError(f"Could not read configuration file: {args.config}")

    PATH_METAFLYE = config["PATHS"]["metaflye"]

    PATH_METAMDBG = config["PATHS"]["metaMDBG"]

    PATH_MINIMAP = config["PATHS"]["minimap"]

    PATH_HAIRSPLITTER = config["PATHS"]["hairsplitter"]

    PATH_MYLOASM = None

    if args.use_myloasm:

        PATH_MYLOASM = config["PATHS"]["myloasm"]

    PREDEFINED_DB = config["DATABASE"]["predefined_db"]

    if config.has_option("DATABASE", "strain_species_json"):

        logging.warning("config.ini: strain_species_json is no longer used (species come from --nodes); ignored.")

    # -------------------------------------------------------------------------
    # Input checks
    # -------------------------------------------------------------------------

    require_file(args.reads, "Reads file")

    require_file(PREDEFINED_DB, "Reference database")

    require_file(args.nodes, "Taxonomy nodes.dmp")

    require_file(args.names, "Taxonomy names.dmp")

    os.makedirs(args.out_folder, exist_ok=True)

    logging.info("Output folder: %s", args.out_folder)

    logging.info("Reads: %s", args.reads)

    logging.info("Read technology: %s", args.reads_flag)

    logging.info("Database: %s", PREDEFINED_DB)

    PREDEFINED_DB_JSON = write_strain_species_json(args.reference_genome_map, args.nodes, os.path.join(args.out_folder, "strain_species.json"))

    logging.info("Threads: %d", args.threads)

    reads_tech_minimap = {"pacbio": ["-c", "-x", "map-pb"], "hifi": ["-c", "-x", "map-hifi"], "ont": ["-c", "-x", "map-ont"]}

    # =========================================================================
    # STEP 1 — Assembly
    # =========================================================================

    logging.info("STEP 1: Assembly")

    assembly_result = run_assembly(args, PATH_METAFLYE, PATH_METAMDBG, PATH_MYLOASM)

    assembly_out_dir = assembly_result["output_dir"]

    assembly_fasta = assembly_result["assembly"]

    assembly_info = assembly_result["assembly_info"]

    assembly_gfa = assembly_result["gfa"]

    assembler = assembly_result["assembler"]

    contig_paths = assembly_result["contig_paths"]

    # =========================================================================
    # STEP 2 — Map assembly to full database
    # =========================================================================

    logging.info("STEP 2: Mapping assembly to full reference database")

    assembly_paf = os.path.join(args.out_folder, "assembly.to_big_db.paf")

    map_assembly_to_database(
        PATH_MINIMAP,
        PREDEFINED_DB,
        assembly_fasta,
        assembly_paf,
        args.threads,
        args.db_map_jobs,
        force=args.force,
        stderr_path=os.path.join(args.out_folder, "minimap_reduction.log"),
    )

    require_file(assembly_paf, "Assembly-to-reference PAF")

    # =========================================================================
    # STEP 3 — HairSplitter
    # =========================================================================

    logging.info("STEP 3: Detecting collapsed strains with HairSplitter")

    hairsplitter_out_dir = os.path.join(assembly_out_dir, "hairsplitter")

    cleaned_assembly = os.path.join(hairsplitter_out_dir, "tmp", "cleaned_assembly.fasta")

    cmd = [PATH_HAIRSPLITTER, "-i", assembly_fasta, "-f", args.reads, "-t", str(args.threads), "-o", hairsplitter_out_dir, "-F"]

    run_command(cmd, target_file=cleaned_assembly, force=args.force, stderr_path=os.path.join(args.out_folder, "hairsplitter.log"))

    collapsed_strains = os.path.join(assembly_out_dir, "collapsed_strains.txt")

    if not os.path.exists(collapsed_strains) or args.force:

        logging.info("Creating collapsed_strains.txt")

        with open(cleaned_assembly) as inp, open(collapsed_strains, "w") as out:

            for line in inp:

                if not line.startswith(">"):
                    continue

                contig = line[1:].strip()

                contig = contig.replace("@", ":")

                contig = contig.replace("edge", "contig")

                out.write(contig + "\n")

    require_file(collapsed_strains, "Collapsed strain information")

    # =========================================================================
    # STEP 4 — Database reduction + novelty classification
    # =========================================================================

    logging.info("STEP 4: Database reduction and contig novelty classification")

    dr_args = Namespace(
        assembly=assembly_fasta,
        contig_paf=assembly_paf,
        strain_species_json=PREDEFINED_DB_JSON,
        reference_genome_map=args.reference_genome_map,
        num_collapsed_strains=collapsed_strains,
        output_dir=args.out_folder,
        database=None,  # the reduced database itself is extracted in step 6, only when the extended database needs (re)building
        threads=args.threads,
        strictness=args.strictness,
        min_contig_len=args.min_contig_len,
        collapsed_strains_overhead=args.collapsed_strains_overhead,
        known_thr=args.known_thr,
        novel_thr=args.novel_thr,
        novelty_score=args.novelty_score,
        no_mapq_filter=args.no_mapq_filter,
    )

    dr_result = DR.run(dr_args)

    classification = dr_result["classification"]

    contig_scores = dr_result["contig_scores"]

    contig_lengths = dr_result["contig_lengths"]

    reduced_refs_path = dr_result["reduced_refs_path"]

    extra_sequences = dr_result["extra_sequences"]

    require_file(reduced_refs_path, "Reduced reference list")

    require_file(extra_sequences, "Extra uncertain/novel contigs")

    # =========================================================================
    # STEP 5 — Novel grouping
    # =========================================================================

    logging.info("STEP 5: Grouping uncertain and novel candidate contigs")

    novel_args = Namespace(
        assembly=assembly_fasta,
        extra_sequences=extra_sequences,
        classification=classification,
        contig_scores=contig_scores,
        contig_lengths=contig_lengths,
        strain_species_json=PREDEFINED_DB_JSON,
        output_dir=args.out_folder,
        threads=args.threads,
        minimap=PATH_MINIMAP,
        force=args.force,
        assembler=assembler,
        assembly_info=assembly_info,
        gfa=assembly_gfa,
        contig_paths=contig_paths,
        grouping_strict_overlap=args.grouping_strict_overlap,
        grouping_strict_identity=args.grouping_strict_identity,
        grouping_relaxed_overlap=args.grouping_relaxed_overlap,
        grouping_relaxed_identity=args.grouping_relaxed_identity,
    )

    novel_result = NovelDetection.run(novel_args)

    extended_species_json = novel_result["extended_strain_species"]

    require_file(extended_species_json, "Extended strain-species JSON")

    # =========================================================================
    # STEP 6 — Build extended database
    # =========================================================================

    logging.info("STEP 6: Building extended database")

    extended_db = os.path.join(args.out_folder, "extended_db.fa")

    host_fa = None

    if args.host_reference:
        host_fa, host_group = prepare_host_reference(args.host_reference, args.host_name, args.out_folder, force=args.force)

        with open(extended_species_json) as fh:
            species_info = json.load(fh)
        if species_info.get(host_group) != host_group:
            species_info[host_group] = host_group
            with open(extended_species_json, "w") as fh:
                json.dump(species_info, fh)
            logging.info("Registered host group %s in %s", host_group, extended_species_json)

    # extended_db.fa = reduced database (reduced_db_refs.txt extracted from the
    # full database) + UNCERTAIN/novel contigs (+ host). A record of what it was
    # built from is kept next to it; when it matches, the extended database and
    # the read mapping against it are reused, otherwise both are rebuilt.
    inputs_record = extended_db + ".inputs"
    host_info = f"{os.path.abspath(host_fa)} {os.path.getsize(host_fa)}" if host_fa else "none"
    inputs_text = (
        f"database\t{os.path.abspath(PREDEFINED_DB)}\n"
        f"reduced_refs_sha256\t{file_sha256(reduced_refs_path)}\n"
        f"extra_sequences_sha256\t{file_sha256(extra_sequences)}\n"
        f"host\t{host_info}\n"
    )
    previous = open(inputs_record).read() if os.path.exists(inputs_record) else None
    rebuild_extended = args.force or not os.path.exists(extended_db) or previous != inputs_text

    if rebuild_extended:
        if os.path.exists(extended_db) and not args.force:
            logging.info("Extended database inputs changed; rebuilding it and the read mapping.")
        reduced_db = os.path.join(args.out_folder, "reduced_db.fa")
        DR.build_reduced_database(PREDEFINED_DB, reduced_refs_path, reduced_db, args.threads)
        concatenate_fastas([reduced_db, extra_sequences] + ([host_fa] if host_fa else []), extended_db, force=True)
        if args.keep_reduced_db:
            logging.info("Reduced database kept: %s", reduced_db)
        else:
            os.remove(reduced_db)
            logging.info(
                "Reduced database removed (its sequences are in extended_db.fa; recreate with: seqkit grep -f %s %s).", reduced_refs_path, PREDEFINED_DB
            )
        with open(inputs_record, "w") as fh:
            fh.write(inputs_text)
    else:
        logging.info("Extended database is up to date with its inputs; reusing it: %s", extended_db)

    require_file(extended_db, "Extended reference database")

    # =========================================================================
    # STEP 7 — Map reads to extended database
    # =========================================================================

    logging.info("STEP 7: Mapping reads to extended database")

    reads_extended_paf = os.path.join(args.out_folder, "reads.to_extended.paf")

    cmd = [PATH_MINIMAP] + reads_tech_minimap[args.reads_flag] + [extended_db, args.reads, "-t", str(args.threads)]

    run_command_stdout(cmd, reads_extended_paf, force=rebuild_extended, stderr_path=os.path.join(args.out_folder, "minimap_classification.log"))

    require_file(reads_extended_paf, "Read-to-extended-database PAF")

    # =========================================================================
    # STEP 8 — StrainRefine
    # =========================================================================

    logging.info("STEP 8: Running StrainRefine")

    read_class_output = os.path.join(args.out_folder, "read_classification.out")

    clustering_out = os.path.join(args.out_folder, "strainrefine_clustering")

    ref_summary_output = os.path.join(args.out_folder, "strainrefine_reference_summary.tsv")

    sr_args = Namespace(
        paf_path=reads_extended_paf,
        strain_species_info=extended_species_json,
        cluster_eps=args.cluster_eps,
        high_score_threshold=args.high_score_threshold,
        min_high_score_reads=args.min_high_score_reads,
        mid_high_score_reads=args.mid_high_score_reads,
        min_high_score_fraction=args.min_high_score_fraction,
        species_min_read_count=args.species_min_read_count,
        species_min_mean_score=args.species_min_mean_score,
        species_low_count_cap=args.species_low_count_cap,
        read_class_output=read_class_output,
        clustering_out=clustering_out,
        ref_summary_output=ref_summary_output,
    )

    sr_args.reference_fasta = extended_db
    sr_args.reference_genome_map = args.reference_genome_map
    sr_args.taxonomy_nodes = args.nodes
    sr_args.taxonomy_names = args.names
    sr_args.threads = args.threads
    sr_args.novel_placement = not args.no_novel_placement
    sr_args.placement_report = args.placement_report
    sr_args.genome_depth = args.genome_depth
    sr_args.min_abundance = args.min_abundance
    sr_args.contig_classification = os.path.join(args.out_folder, "contig_classification.tsv")
    sr_args.taxonomic_profile_output = None

    sr_result = StrainRefine.run(sr_args)

    require_file(read_class_output, "StrainRefine read classification")

    logging.info("StrainRefine finished successfully.")

    # =========================================================================
    # STEP 9 — Abundance
    # =========================================================================

    logging.info("STEP 9: Calculating abundances")

    CalculateAbundances.run(sr_result, sr_args)

    require_file(ref_summary_output, "Reference summary")

    logging.info("Abundance calculation finished successfully.")

    logging.info("MADRe finished successfully.")


if __name__ == "__main__":
    main()
