#!/usr/bin/env python3

"""
NovelDetection
=================

Novel-organism grouping stage for MADRe (step 5 of MADRe).
Taxonomic labels for novel contigs are not given here: StrainRefine
places each novel contig in the taxonomy with read evidence.

Input is produced by DatabaseReduction:

    - extra_sequences.fa containing UNCERTAIN + NOVEL_CANDIDATE contigs
    - contig classification
    - contig scores / lengths
    - strain -> species JSON

Grouping:

    1. all-vs-all mapping of extra contigs
    2. strict sequence-similarity groups
    3. relaxed sequence-similarity edges
    4. assembly-graph edges (metaFlye, Myloasm, metaMDBG) and Myloasm
       alternate -> primary contig edges
    5. connected comparison groups

Final StrainRefine mapping:

    KNOWN database reference
        -> its species from the strain -> species JSON

    UNCERTAIN assembly contig
        -> species of its best known reference

    NOVEL_CANDIDATE assembly contig
        -> NOVEL_COMPARE_XXXXX synthetic species
"""

import json
import logging
import os
import re
import subprocess
from collections import Counter, defaultdict


# =============================================================================
# Defaults
# =============================================================================

NOVEL_PREFIX = "NOVEL"
NOVEL_COMPARE_PREFIX = "NOVEL_COMPARE"


# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[logging.StreamHandler()])


# =============================================================================
# Helpers
# =============================================================================


def norm_id(x):
    """
    Normalise assembly contig IDs ("scaffold" -> "contig").
    """
    return x.replace("scaffold", "contig")


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def load_fasta_ids(fasta):
    ids = set()

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(norm_id(line[1:].split()[0]))

    return ids


def extract_taxid(reference):
    """
    Extract taxid from a MADRe reference ID such as:

        |12345|NZ_CP000001.1

    or:

        taxid|12345|NZ_CP000001.1
    """

    if not reference:
        return None

    m = re.search(r"(?:taxid)?\|(\d+)\|", reference)

    if m:
        return m.group(1)

    parts = reference.split("|")

    if len(parts) > 1 and parts[1].isdigit():
        return parts[1]

    return None


# =============================================================================
# Union-find
# =============================================================================


class UnionFind:

    def __init__(self, nodes):
        self.parent = {n: n for n in nodes}

        self.rank = {n: 0 for n in nodes}

    def find(self, x):

        while self.parent[x] != x:

            self.parent[x] = self.parent[self.parent[x]]

            x = self.parent[x]

        return x

    def union(self, a, b):

        ra = self.find(a)
        rb = self.find(b)

        if ra == rb:
            return False

        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra

        self.parent[rb] = ra

        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

        return True

    def components(self):

        comp = defaultdict(set)

        for n in self.parent:
            comp[self.find(n)].add(n)

        return list(comp.values())


# =============================================================================
# Assembly information: coverage, circularity and contig links
#
# Every assembler gives the same three things:
#     coverage      contig -> coverage
#     circular      contig -> "yes" / "no"
#     edge_contigs  graph segment -> set of contigs made of that segment,
#                   so the GFA links between segments become contig links
#
#     metaFlye  assembly_info.txt (coverage, circular, contig -> graph edges) + assembly_graph.gfa
#     Myloasm   FASTA headers (depth, circular, alternate -> primary contig) + final_contig_graph.gfa
#               (its segments are the contigs, named by the header part before "_len-")
#     metaMDBG  FASTA headers (coverage, circular) + "metaMDBG gfa" output: contig -> unitig table
#               (assemblyGraph_k*_contigPath.tsv) + graph without sequences (assemblyGraph_k*.noseq.gfa)
# =============================================================================


def read_fasta_headers(fasta):
    """contig -> rest of its FASTA header line (after the name)."""
    headers = {}
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name, _, rest = line[1:].rstrip("\n").partition(" ")
                headers[norm_id(name)] = rest
    return headers


def load_flye_info(assembly_info):
    """Flye's assembly_info.txt: coverage, circular and contig -> graph edges ("edge_N")."""
    cov, circular, edge_contigs = {}, {}, defaultdict(set)
    with open(assembly_info) as fh:
        next(fh, None)
        for line in fh:
            p = line.strip().split("\t")
            if len(p) < 8:
                continue
            cid = norm_id(p[0])
            try:
                cov[cid] = float(p[2])
            except ValueError:
                cov[cid] = 0.0
            circular[cid] = "yes" if p[3] == "Y" else "no"
            if p[7] in {"*", ""}:
                continue
            for token in p[7].split(","):
                edge = token.lstrip("-")
                if edge not in {"*", ""}:
                    edge_contigs["edge_" + edge].add(cid)
    return cov, circular, edge_contigs


def load_myloasm_info(assembly_fasta):
    """Myloasm FASTA headers: depth (first value), circular, graph segment name, and alternate -> primary links."""
    cov, circular, edge_contigs, alternate_links = {}, {}, defaultdict(set), []
    names = read_fasta_headers(assembly_fasta)
    for cid, rest in names.items():
        m = re.search(r"_depth-([0-9.]+)", cid)
        if m:
            cov[cid] = float(m.group(1))
        m = re.search(r"_circular-(yes|no)", cid)
        if m:
            circular[cid] = m.group(1)
        edge_contigs[cid.split("_len-")[0]].add(cid)
        m = re.search(r"secondary_to:\|(\S+)", rest)
        if m and norm_id(m.group(1)) in names:
            alternate_links.append((cid, norm_id(m.group(1))))
    return cov, circular, edge_contigs, alternate_links


def load_metamdbg_info(assembly_fasta, contig_paths):
    """metaMDBG FASTA headers (coverage=, circular=) and the contig -> unitig table of "metaMDBG gfa"."""
    cov, circular, edge_contigs = {}, {}, defaultdict(set)
    for cid, rest in read_fasta_headers(assembly_fasta).items():
        m = re.search(r"coverage=([0-9.]+)", rest)
        if m:
            cov[cid] = float(m.group(1))
        m = re.search(r"circular=(yes|no)", rest)
        if m:
            circular[cid] = m.group(1)
    if contig_paths and os.path.exists(contig_paths):
        with open(contig_paths) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    edge_contigs[p[1]].add(norm_id(p[0]))
    return cov, circular, edge_contigs


def load_assembly_info(assembler, assembly_fasta, assembly_info=None, contig_paths=None):
    """
    Coverage, circularity, graph segment -> contigs, and alternate -> primary contig links
    for the assembler that was used. Missing files give empty results.
    """
    if assembler == "metaflye":
        if assembly_info and os.path.exists(assembly_info):
            return load_flye_info(assembly_info) + ([],)
    elif assembler == "myloasm":
        return load_myloasm_info(assembly_fasta)
    elif assembler == "metamdbg":
        return load_metamdbg_info(assembly_fasta, contig_paths) + ([],)
    elif assembler is not None:
        raise ValueError(f"Unknown assembler: {assembler}")
    return {}, {}, defaultdict(set), []


# =============================================================================
# Assembly graph
# =============================================================================


def load_gfa_edges(gfa, edge_contigs, valid_contigs):
    """
    Convert GFA links between graph segments into links between assembly contigs.
    """
    graph_edges = []
    if gfa is None or not os.path.exists(gfa):
        return graph_edges
    with open(gfa) as fh:
        for line in fh:
            if not line.startswith("L"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            for contig_a in edge_contigs.get(p[1], set()):
                for contig_b in edge_contigs.get(p[3], set()):
                    if contig_a != contig_b and contig_a in valid_contigs and contig_b in valid_contigs:
                        graph_edges.append((contig_a, contig_b))
    return graph_edges


# =============================================================================
# All-vs-all contig mapping
# =============================================================================


def run_contig_self_alignment(minimap_path, extra_fasta, output_paf, threads, force=False):
    """
    Run minimap2 all-vs-all alignment of extra assembly contigs.
    """

    if os.path.exists(output_paf) and not force:

        logging.info("Using existing extra-contig PAF: %s", output_paf)

        return

    cmd = [minimap_path, "-x", "asm5", "-t", str(threads), extra_fasta, extra_fasta]

    logging.info("Running extra-contig all-vs-all mapping: %s", " ".join(cmd))

    with open(output_paf, "w") as out:

        subprocess.run(cmd, stdout=out, check=True)


# =============================================================================
# PAF grouping edges
# =============================================================================


def parse_contig_paf_edges(paf, valid_contigs, strict_overlap, strict_identity, relaxed_overlap, relaxed_identity):
    """
    Overlap and identity of each contig-contig alignment:

        query overlap  = aligned query span / query length
        target overlap = aligned target span / target length

        overlap = max(query overlap, target overlap)

        identity = matches / PAF alignment block length
    """

    strict_edges = []
    relaxed_edges = []

    with open(paf) as fh:

        for line in fh:

            p = line.rstrip("\n").split("\t")

            if len(p) < 12:
                continue

            qid = norm_id(p[0])

            tid = norm_id(p[5])

            if qid == tid:
                continue

            if qid not in valid_contigs or tid not in valid_contigs:
                continue

            qlen, qstart, qend = int(p[1]), int(p[2]), int(p[3])
            tlen, tstart, tend, matches, block_len = int(p[6]), int(p[7]), int(p[8]), int(p[9]), int(p[10])

            qov = (qend - qstart) / qlen if qlen else 0.0

            tov = (tend - tstart) / tlen if tlen else 0.0

            overlap = max(qov, tov)

            identity = matches / block_len if block_len else 0.0

            if overlap >= strict_overlap and identity >= strict_identity:

                strict_edges.append((qid, tid))

            if overlap >= relaxed_overlap and identity >= relaxed_identity:

                relaxed_edges.append((qid, tid))

    return (strict_edges, relaxed_edges)


# =============================================================================
# Group IDs
# =============================================================================


def groups_to_map(groups, prefix):
    result = {}

    groups = sorted(groups, key=lambda g: sorted(g)[0])

    for i, group in enumerate(groups, start=1):

        gid = f"{prefix}_{i:05d}"

        for cid in group:
            result[cid] = gid

    return result


# =============================================================================
# UNCERTAIN assignments
# =============================================================================


def assign_uncertain_species(classification, strain_species):
    """
    UNCERTAIN contigs are assigned to the species corresponding to their
    best known database reference.
    """

    assignments = {}

    for cid, (label, _score, best_ref) in classification.items():

        if label != "UNCERTAIN":
            continue

        taxid = extract_taxid(best_ref)

        if taxid is None:

            logging.warning("Could not extract taxid from best reference for " "UNCERTAIN contig %s: %s", cid, best_ref)

            continue

        species_taxid = strain_species.get(taxid)

        if species_taxid is None:

            logging.warning("Taxid %s for UNCERTAIN contig %s is not present " "in strain-species JSON.", taxid, cid)

            continue

        assignments[norm_id(cid)] = species_taxid

    return assignments


# =============================================================================
# Extended StrainRefine taxonomy mapping
# =============================================================================


def build_extended_json(original_json, classification, comparison_map, output_json):
    """
    Construct the mapping used by StrainRefine.

    Original database references:
        unchanged

    UNCERTAIN contigs:
        assigned to species of best known reference

    NOVEL_CANDIDATE contigs:
        assigned to synthetic NOVEL_COMPARE group
    """

    strain_species = load_json(original_json)

    uncertain_species = assign_uncertain_species(classification, strain_species)

    extended = dict(strain_species)

    # -------------------------------------------------------------------------
    # UNCERTAIN -> known species
    # -------------------------------------------------------------------------

    for cid, species_taxid in uncertain_species.items():

        extended[cid] = species_taxid

    # -------------------------------------------------------------------------
    # NOVEL_CANDIDATE -> synthetic comparison species
    # -------------------------------------------------------------------------

    missing_novel = []

    for cid, (label, _score, _best_ref) in classification.items():

        cid = norm_id(cid)

        if label != "NOVEL_CANDIDATE":
            continue

        if cid not in comparison_map:

            missing_novel.append(cid)

            continue

        extended[cid] = comparison_map[cid]

    if missing_novel:

        logging.warning("%d NOVEL_CANDIDATE contigs have no comparison group.", len(missing_novel))

    write_json(output_json, extended)

    return (extended, uncertain_species)


# =============================================================================
# Reports
# =============================================================================


def write_group_summary(out_tsv, contigs, classification, strict_map, comparison_map, method, coverage, circular):
    """
    Integrated grouping report.
    """

    with open(out_tsv, "w") as fh:

        fh.write("contig_id\t" "classification\t" "vcig_norm\t" "strict_novel_group\t" "comparison_group\t" "method\t" "coverage\t" "circular\t" "best_ref\n")

        for cid in sorted(contigs):

            class_values = classification.get(cid, classification.get(norm_id(cid), ("", 0.0, None)))

            label, vcig, best_ref = class_values

            fh.write(
                f"{cid}\t"
                f"{label}\t"
                f"{vcig:.6f}\t"
                f"{strict_map.get(cid, '')}\t"
                f"{comparison_map.get(cid, '')}\t"
                f"{method.get(cid, '')}\t"
                f"{coverage.get(cid, '')}\t"
                f"{circular.get(cid, '')}\t"
                f"{best_ref or ''}\n"
            )


# =============================================================================
# Main NovelDetection API
# =============================================================================


def run(args):
    """
    Run MADRe novel grouping.

    Expected args
    -------------
    extra_sequences
    classification
    contig_scores
    contig_lengths

    strain_species_json

    output_dir
    threads

    minimap
    force

    assembly
    assembler         metaflye / myloasm / metamdbg (None: no assembly information)
    assembly_info     metaFlye assembly_info.txt
    gfa               assembly graph (metaFlye, Myloasm, or metaMDBG without sequences)
    contig_paths      metaMDBG contig -> unitig table

    grouping_strict_overlap
    grouping_strict_identity
    grouping_relaxed_overlap
    grouping_relaxed_identity

    Returns
    -------
    dict
        extended_strain_species   path of extended_strain_species.json (all other outputs are files in output_dir)
    """

    os.makedirs(args.output_dir, exist_ok=True)

    logging.info("Starting NovelDetection")

    # =========================================================================
    # STEP 1 — Load extra contigs
    # =========================================================================

    contigs = load_fasta_ids(args.extra_sequences)

    logging.info("Extra contigs: %d", len(contigs))

    if not contigs:

        logging.warning("No UNCERTAIN or NOVEL_CANDIDATE contigs detected.")

    # =========================================================================
    # STEP 2 — Assembly coverage and graph information
    # =========================================================================

    assembler = getattr(args, "assembler", None)
    coverage, circular, edge_contigs, alternate_links = load_assembly_info(
        assembler, args.assembly, getattr(args, "assembly_info", None), getattr(args, "contig_paths", None)
    )
    graph_edges = load_gfa_edges(getattr(args, "gfa", None), edge_contigs, contigs)
    alternate_edges = [(a, b) for a, b in alternate_links if a in contigs and b in contigs]

    logging.info("Assembler: %s | contigs with coverage: %d | circular contigs: %d", assembler, len(coverage), sum(v == "yes" for v in circular.values()))
    logging.info("Assembly graph edges between extra contigs: %d | alternate -> primary edges: %d", len(graph_edges), len(alternate_edges))

    # =========================================================================
    # STEP 3 — Extra contig all-vs-all alignment
    # =========================================================================

    contig_paf = os.path.join(args.output_dir, "extra_sequences.self.paf")

    if contigs:

        run_contig_self_alignment(
            minimap_path=args.minimap, extra_fasta=args.extra_sequences, output_paf=contig_paf, threads=args.threads, force=getattr(args, "force", False)
        )

    else:

        with open(contig_paf, "w"):
            pass

    # =========================================================================
    # STEP 4 — Sequence-similarity edges
    # =========================================================================

    strict_edges, relaxed_edges = parse_contig_paf_edges(
        contig_paf, contigs, args.grouping_strict_overlap, args.grouping_strict_identity, args.grouping_relaxed_overlap, args.grouping_relaxed_identity
    )

    logging.info("Strict PAF edges: %d", len(strict_edges))

    logging.info("Relaxed PAF edges: %d", len(relaxed_edges))

    # =========================================================================
    # STEP 5 — Strict NOVEL groups
    # =========================================================================

    uf_strict = UnionFind(contigs)

    for a, b in strict_edges:

        uf_strict.union(a, b)

    strict_map = groups_to_map(uf_strict.components(), NOVEL_PREFIX)

    # =========================================================================
    # STEP 6 — NOVEL_COMPARE groups: joined by strict, relaxed, assembly-graph and alternate -> primary edges
    # =========================================================================

    uf_compare = UnionFind(contigs)

    for a, b in strict_edges + relaxed_edges + graph_edges + alternate_edges:

        uf_compare.union(a, b)

    comparison_map = groups_to_map(uf_compare.components(), NOVEL_COMPARE_PREFIX)

    # =========================================================================
    # STEP 7 — Group membership method (for novel_groups.tsv): "connected"
    # if the contig shares its comparison group with other contigs,
    # "singleton" if it is alone.
    # =========================================================================

    group_size = defaultdict(int)
    for gid in comparison_map.values():
        group_size[gid] += 1

    method = {cid: ("connected" if group_size[comparison_map[cid]] > 1 else "singleton") for cid in contigs}

    # =========================================================================
    # STEP 8 — Output grouping maps
    # =========================================================================

    strict_json = os.path.join(args.output_dir, "novel_id_map.json")

    comparison_json = os.path.join(args.output_dir, "novel_comparison_group_map.json")

    summary_tsv = os.path.join(args.output_dir, "novel_groups.tsv")

    write_json(strict_json, strict_map)

    write_json(comparison_json, comparison_map)

    write_group_summary(summary_tsv, contigs, args.classification, strict_map, comparison_map, method, coverage, circular)

    # =========================================================================
    # STEP 9 — Extended strain/species mapping
    # =========================================================================

    extended_json = os.path.join(args.output_dir, "extended_strain_species.json")

    (extended_mapping, uncertain_species) = build_extended_json(args.strain_species_json, args.classification, comparison_map, extended_json)

    # =========================================================================
    # Summary
    # =========================================================================

    n_uncertain = sum(1 for label, _, _ in args.classification.values() if label == "UNCERTAIN")

    n_novel = sum(1 for label, _, _ in args.classification.values() if label == "NOVEL_CANDIDATE")

    logging.info("NovelDetection complete.")

    logging.info("UNCERTAIN contigs: %d", n_uncertain)

    logging.info("NOVEL_CANDIDATE contigs: %d", n_novel)

    logging.info("Strict groups: %d", len(set(strict_map.values())))

    logging.info("Comparison groups: %d", len(set(comparison_map.values())))

    logging.info("UNCERTAIN contigs assigned to known species: %d", len(uncertain_species))

    return {"extended_strain_species": extended_json}
