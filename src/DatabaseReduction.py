#!/usr/bin/env python3
"""
DatabaseReduction
========================

Database reduction for MADRe (step 4 of MADRe), with genome-atomic
reduced-reference selection:

    contig mappings
        -> U / NU classification
        -> NU cleaning
        -> EM
        -> species-aware reduced reference selection
        -> genome-atomic expansion

Genome-atomic expansion: once any sequence belonging to a source genome
(genome_id, e.g. a GCF/GCA assembly accession) is selected for the reduced
database, every other sequence of that genome is added too, using
reference_genome_map.tsv to look up which sequences belong to which genome.
Without it, a multi-contig draft genome could enter the reduced database only
partially (e.g. 11 of its 58 contigs), letting a more complete competing
reference win the reads that truly belong to it. The selection decisions
(U/NU/EM/species logic) do not depend on it; only the final set of references
written to the reduced database is expanded to whole genomes.

MADRe additionally classifies assembly contigs using normalized vcig
scores:

    KNOWN
    UNCERTAIN
    NOVEL_CANDIDATE

UNCERTAIN and NOVEL_CANDIDATE contigs are extracted to extra_sequences.fa
for later novel grouping and construction of the extended reference
database. This classification and extraction are unaffected by the
genome-atomic expansion, which only applies to the reduced (known-species)
reference set.

The module does NOT perform:
    - novel contig grouping
    - extended database construction
    - read mapping
    - StrainRefine

Those steps belong to later stages of MADRe.
"""

import argparse
import json
import logging
import os
import subprocess
from collections import defaultdict


# =============================================================================
# Constants
# =============================================================================

CLEAN_THRESHOLD = {"less-strict": 0.9, "strict": 0.99, "very-strict": 0.999}

MIN_CONTIG_LEN = 1000
MAX_COLLAPSED_STRAINS_OVERHEAD = 2

DEFAULT_KNOWN_THR = 0.75
DEFAULT_NOVEL_THR = 0.35

EM_MAX_ITER = 25
EM_EPSILON = 0.0001


# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[logging.StreamHandler()])


# =============================================================================
# General helpers
# =============================================================================


def norm_contig_id(contig_id):
    """
    Normalise assembly contig IDs ("scaffold" -> "contig").
    """
    return contig_id.replace("scaffold", "contig")


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_reference_list(references, output_path):
    references = sorted(set(references))

    with open(output_path, "w") as fh:
        for ref in references:
            fh.write(ref + "\n")

    logging.info("Reduced reference list: %d references -> %s", len(references), output_path)


# =============================================================================
# Reference genome map (contig -> source genome)
# =============================================================================


def load_reference_genome_map(path):
    """
    Load reference_genome_map.tsv (header: reference, genome_id,
    genome_length) built by PrepareDatabase.py.

    Returns
    -------
    ref_to_genome
        Reference header -> genome_id (e.g. a GCF/GCA assembly accession).

    genome_to_refs
        genome_id -> set of every reference header belonging to it.
    """

    ref_to_genome = {}
    genome_to_refs = defaultdict(set)

    with open(path) as fh:

        header = fh.readline()

        if not header:
            raise ValueError(f"Empty reference genome map: {path}")

        for line in fh:

            line = line.rstrip("\n")

            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < 2:
                logging.warning("Skipping malformed reference genome map line: %s", line)
                continue

            reference, genome_id = parts[0], parts[1]

            ref_to_genome[reference] = genome_id
            genome_to_refs[genome_id].add(reference)

    logging.info("Loaded reference genome map: %d references across %d genomes -> %s", len(ref_to_genome), len(genome_to_refs), path)

    return ref_to_genome, genome_to_refs


def expand_to_genome_level(references, ref_to_genome, genome_to_refs):
    """
    Genome-atomic expansion of a selected reference set.

    If any reference belonging to a genome_id is present in `references`,
    every reference belonging to that same genome_id is added too.
    References absent from the reference genome map are kept as-is,
    unexpanded, with a warning.
    """

    expanded = set(references)

    unmapped = [ref for ref in references if ref not in ref_to_genome]

    if unmapped:

        logging.warning(
            "%d selected reference(s) not present in the reference genome " "map; kept without genome-level expansion (showing up to 5): %s",
            len(unmapped),
            unmapped[:5],
        )

    genome_ids = {ref_to_genome[ref] for ref in references if ref in ref_to_genome}

    added = 0

    for genome_id in genome_ids:

        for ref in genome_to_refs.get(genome_id, ()):

            if ref not in expanded:

                expanded.add(ref)
                added += 1

    if added:

        logging.info(
            "Genome-level expansion: added %d sibling reference(s) across " "%d genome(s) so no selected genome is left partially " "represented.",
            added,
            len(genome_ids),
        )

    else:

        logging.info("Genome-level expansion: all selected genomes were already " "fully represented; nothing added.")

    return sorted(expanded)


# =============================================================================
# MADRe mapping classification
# =============================================================================


def _mapping_output(mapping_class_path, mapping_reduced_path, predictions_mapping):
    """
    Write MADRe-compatible contig mapping outputs.
    """

    reduced = []

    if mapping_class_path is not None:
        with open(mapping_class_path, "w") as fh:
            for contig, ref in predictions_mapping.items():
                reduced.append(ref)
                fh.write(f"{contig} : {ref}\n")

    if mapping_reduced_path is not None:

        if mapping_class_path is None:
            reduced.extend(predictions_mapping.values())

        with open(mapping_reduced_path, "w") as fh:
            for ref in sorted(set(reduced)):
                fh.write(ref + "\n")

    logging.info("Mapping information written.")


def read_contig_paf(paf_path, min_contig_len=MIN_CONTIG_LEN, no_mapq_filter=False):
    """
    Read the contig-to-reference PAF once.

    Skips incomplete lines, contigs shorter than min_contig_len, MAPQ=0
    alignments (unless no_mapq_filter) and alignments without matches. For
    every remaining alignment:

        vcig = 2 * length * matches / (length + matches)
        length = max(aligned query span, aligned target span)

    Returns
    -------
    vcig
        contig -> {reference: vcig summed over the contig's alignments to it}
        (contigs and references in order of first appearance in the PAF)
    contig_lengths
        contig -> length, for every contig in the PAF (also short ones)
    references
        references in order of first appearance among the kept alignments
    alignments
        contig -> {reference: [(query start, query end, vcig per aligned base)]},
        used by the "union" novelty score
    """
    vcig = {}
    contig_lengths = {}
    references = {}
    alignments = {}
    n_lines = n_malformed = n_short = n_mapq0 = n_zero = 0

    with open(paf_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 12:
                n_malformed += 1
                continue
            n_lines += 1
            contig = norm_contig_id(parts[0])
            qlen = int(parts[1])
            contig_lengths[contig] = qlen
            if qlen < min_contig_len:
                n_short += 1
                continue
            if int(parts[11]) == 0 and not no_mapq_filter:
                n_mapq0 += 1
                continue
            length = max(int(parts[3]) - int(parts[2]), int(parts[8]) - int(parts[7]))
            matches = int(parts[9])
            if matches <= 0:
                n_zero += 1
                continue
            if length <= 0:
                continue
            ref = parts[5]
            references.setdefault(ref, len(references))
            per_contig = vcig.setdefault(contig, {})
            per_contig[ref] = per_contig.get(ref, 0.0) + (2.0 * length * matches) / (length + matches)
            alignments.setdefault(contig, {}).setdefault(ref, []).append((int(parts[2]), int(parts[3]), 2.0 * matches / (length + matches)))

    if n_malformed:
        logging.warning(f"Skipped {n_malformed} malformed PAF lines with fewer than 12 fields.")
    logging.info(f"Contig PAF: {n_lines} alignments; {n_short} from short contigs; {n_mapq0} MAPQ=0 filtered; {n_zero} zero-match.")
    return vcig, contig_lengths, list(references), alignments


def calculate_mapping_class(vcig, references, mapping_class_path=None, mapping_reduced_path=None):
    """
    Construct MADRe U / NU classes from the summed vcig per contig/reference
    pair (read_contig_paf).

    Returns
    -------
    U
        Unique contig mappings.
    NU
        Non-unique contig mappings.
    genomes
        Integer genome ID -> reference identifier mapping.
    """
    logging.info("Extracting mapping information.")
    predictions_mapping = {}
    predictions_mapping_vcg = vcig
    U = {}
    NU = {}
    genomes = dict(enumerate(references))
    genomes_names = {ref: i for i, ref in genomes.items()}

    # -------------------------------------------------------------------------
    # Convert accumulated mappings into U / NU
    # -------------------------------------------------------------------------

    for read_id, candidates in predictions_mapping_vcg.items():

        max_value = max(candidates.values())

        for candidate, value_cig in candidates.items():

            ref_id = genomes_names[candidate]

            if read_id not in U and read_id not in NU:

                U[read_id] = [[ref_id], [value_cig], [float(value_cig)], value_cig]

                predictions_mapping[read_id] = candidate
                continue

            if value_cig == max_value:
                predictions_mapping[read_id] = candidate

            if read_id in U:

                if ref_id in U[read_id][0]:
                    continue

                NU[read_id] = U[read_id]
                del U[read_id]

            if ref_id in NU[read_id][0]:
                continue

            NU[read_id][0].append(ref_id)
            NU[read_id][1].append(value_cig)

            if value_cig > NU[read_id][3]:
                NU[read_id][3] = float(value_cig)

    logging.info("Mapping classification finished: %d unique, %d non-unique contigs.", len(U), len(NU))

    if mapping_class_path is not None or mapping_reduced_path is not None:
        _mapping_output(mapping_class_path, mapping_reduced_path, predictions_mapping)

    return U, NU, genomes


# =============================================================================
# Collapsed-strain information
# =============================================================================


def get_number_of_clusters(num_collapsed_strains):
    """
    Read HairSplitter collapsed-strain information.

    Input format:
        contig_id:number

    MADRe uses number + 1 as the total number of expected clusters.
    """

    number_of_clusters = {}

    with open(num_collapsed_strains) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            parts = line.split(":")

            if len(parts) < 2:
                logging.warning("Skipping malformed collapsed-strain line: %s", line)
                continue

            read_id = norm_contig_id(parts[0])

            number_of_clusters[read_id] = int(parts[1]) + 1

    return number_of_clusters


# =============================================================================
# NU cleaning
# =============================================================================


def clear_NU(U, NU, number_of_clusters, strictness, collapsed_strains_overhead):
    """
    Clean non-unique mappings using MADRe's strictness threshold and
    HairSplitter collapsed-strain estimate.
    """

    NU_cleaned = {}

    for contig, nu_values in NU.items():

        candidate_scores = list(zip(nu_values[0], nu_values[1]))

        max_score = max(nu_values[1])

        # If HairSplitter has no explicit value for the contig, assume one
        # cluster rather than crashing.
        cluster_count = number_of_clusters.get(contig, 1)

        max_num_clusters = max(cluster_count + collapsed_strains_overhead, 1)

        candidate_scores = sorted(candidate_scores, key=lambda x: x[1], reverse=True)[:max_num_clusters]

        cleaned = [(candidate, score) for candidate, score in candidate_scores if score > CLEAN_THRESHOLD[strictness] * max_score]

        # The top candidate should normally always survive, but keep this
        # defensive fallback so the NU structure never becomes empty.
        if not cleaned:
            cleaned = [candidate_scores[0]]

        genome_ids, scores = zip(*cleaned)

        NU_cleaned[contig] = [list(genome_ids), list(scores), nu_values[2], nu_values[3]]

    return U, NU_cleaned


# =============================================================================
# EM
# =============================================================================


def EM(U, NU, genomes, max_iter=EM_MAX_ITER, em_epsilon=EM_EPSILON):
    """
    Original MADRe expectation-maximization procedure.
    """

    logging.info("Starting EM.")

    if len(genomes) == 0:
        logging.warning("No mapped reference genomes available for EM.")
        return [], [], [], NU

    pi = [1 / len(genomes) for _ in genomes]

    init_pi = pi[:]

    pisum_beg = [0 for _ in genomes]

    U_best = [U[i][1][0] for i in U]

    U_total = sum(U_best) if U_best else 0

    for contig in U:

        genome_id = U[contig][0][0]

        pisum_beg[genome_id] += U[contig][1][0] / U[contig][3]

    if U_total > 0:

        pisum_beg = [value / U_total for value in pisum_beg]

    len_NU = 1 if len(NU) == 0 else len(NU)

    thetasum = [0 for _ in genomes]

    iteration = 0

    for iteration in range(max_iter):

        pi_old = pi[:]

        thetasum = [0 for _ in genomes]

        # ---------------------------------------------------------------------
        # E step
        # ---------------------------------------------------------------------

        for contig in NU:

            z = NU[contig]

            indices = z[0]

            pitmp = [pi[k] for k in indices]

            xtmp = [pitmp[k] * z[1][k] / z[3] for k in range(len(indices))]

            xsum = sum(xtmp)

            if xsum == 0:
                xnorm = [0.0 for _ in xtmp]
            else:
                xnorm = [value / xsum for value in xtmp]

            NU[contig][2] = xnorm

            for k in range(len(indices)):
                thetasum[indices[k]] += xnorm[k]

        # ---------------------------------------------------------------------
        # M step
        # ---------------------------------------------------------------------

        pisum = [thetasum[k] + pisum_beg[k] for k in range(len(thetasum))]

        denominator = len(U) + len(NU) + len(pisum)

        if denominator == 0:
            break

        pi = [value / denominator for value in pisum]

        if iteration == 0:
            init_pi = pi[:]

        cutoff = sum(abs(pi_old[k] - pi[k]) for k in range(len(pi)))

        if cutoff <= em_epsilon or len_NU == 1:
            break

    logging.info("Number of performed EM iterations: %d", iteration)

    return init_pi, pi, thetasum, NU


# =============================================================================
# Species-aware reduced-reference selection
# =============================================================================


def get_reduced_list(U, NU, genomes_ids, strain_species_json):
    """
    Construct the reduced reference list using MADRe's species-aware logic.

    This returns per-contig selections only. Genome-atomic expansion is
    applied afterwards, in run(), via expand_to_genome_level().
    """

    logging.info("Finalizing list of references for reduced database.")

    reduced_list = []

    # -------------------------------------------------------------------------
    # Unique mappings
    # -------------------------------------------------------------------------

    for contig, values in U.items():

        ref = genomes_ids[values[0][0]]

        reduced_list.append(ref)

        logging.info("Adding UNIQUE contig %s -> %s", contig, ref)

    # -------------------------------------------------------------------------
    # Non-unique mappings
    # -------------------------------------------------------------------------

    strain_species = load_json(strain_species_json)

    for contig, values in NU.items():

        if not values[2]:
            continue

        max_probability = max(values[2])

        if max_probability == 0:
            logging.warning("All EM probabilities are zero for contig %s. Skipping.", contig)
            continue

        normalized_probabilities = [probability / max_probability for probability in values[2]]

        refs = [genomes_ids[values[0][i]] for i in range(len(values[2]))]

        candidates = sorted(zip(refs, normalized_probabilities), key=lambda x: x[1], reverse=True)

        # ---------------------------------------------------------------------
        # Determine most-supported species
        # ---------------------------------------------------------------------

        species_support = defaultdict(float)

        valid_candidates = []

        for ref, support in candidates:

            parts = ref.split("|")

            if len(parts) < 2:
                logging.warning("Could not extract taxid from reference: %s", ref)
                continue

            taxid = parts[1]

            species_taxid = strain_species.get(taxid)

            if species_taxid is None:

                logging.warning("Taxid %s not present in strain-species JSON. " "Skipping candidate %s for contig %s.", taxid, ref, contig)

                continue

            species_support[species_taxid] += support

            valid_candidates.append((ref, species_taxid, support))

        if not species_support:

            logging.warning("No valid taxonomic candidates for contig %s.", contig)

            continue

        selected_species = max(species_support, key=species_support.get)

        # ---------------------------------------------------------------------
        # Retain all surviving strains belonging to selected species
        # ---------------------------------------------------------------------

        for ref, species_taxid, _ in valid_candidates:

            if species_taxid == selected_species:

                reduced_list.append(ref)

                logging.info("Adding NON-UNIQUE contig %s -> %s", contig, ref)

            else:

                logging.info("Removing contig %s candidate %s", contig, ref)

    return reduced_list


# =============================================================================
# Novel-detection scoring
# =============================================================================


def union_vcig(intervals):
    """Each contig base counted once: the aligned query stretches are merged, and a base is weighted by the
    best alignment covering it (vcig per aligned base, 2M / (L + M))."""
    points = sorted({x for start, end, _ in intervals for x in (start, end)})
    total = 0.0
    for a, b in zip(points, points[1:]):
        weight = max((w for start, end, w in intervals if start <= a and end >= b), default=0.0)
        total += (b - a) * weight
    return total


def contig_novelty_scores(vcig, contig_lengths, alignments=None, method="summed"):
    """
    Normalized vcig scores for each contig/reference pair, from what
    read_contig_paf returns:

        summed: vcig summed over the contig's alignments to the reference / contig length
                (overlapping alignments, e.g. repeats, are all counted, so it can exceed 1)
        union:  each contig base counted once (union_vcig) / contig length, between 0 and 1
    """
    if method not in ("summed", "union"):
        raise ValueError(f"Unknown novelty score: {method}")
    contig_scores = {}
    for cid, ref_vcigs in vcig.items():
        qlen = contig_lengths.get(cid, 0)
        if qlen <= 0:
            continue
        if method == "union":
            contig_scores[cid] = {ref: union_vcig(alignments[cid][ref]) / qlen for ref in ref_vcigs}
        else:
            contig_scores[cid] = {ref: value / qlen for ref, value in ref_vcigs.items()}
    logging.info(f"Novelty scores generated for {len(contig_scores)} contigs.")
    return contig_scores


# =============================================================================
# Contig classification
# =============================================================================


def read_assembly_lengths(assembly_fasta):
    """
    Read contig lengths directly from assembly FASTA.
    """

    lengths = {}

    current_id = None

    with open(assembly_fasta) as fh:

        for line in fh:

            line = line.rstrip()

            if line.startswith(">"):

                current_id = norm_contig_id(line[1:].split()[0])

                lengths[current_id] = 0

            elif current_id is not None:

                lengths[current_id] += len(line)

    return lengths


def classify_contigs(assembly_fasta, contig_scores, paf_contig_lengths, min_contig_len, known_thr, novel_thr, ref_scores=None, selected_refs=None):
    """
    Classify each assembly contig as:

        KNOWN
            best vcig_norm >= known_thr

        UNCERTAIN
            novel_thr <= best vcig_norm < known_thr

        NOVEL_CANDIDATE
            best vcig_norm < novel_thr or no mapping

    The label comes from the highest score in contig_scores. The best reference is the one
    with the highest score, unless ref_scores and selected_refs are given: then it is the
    reference with the highest ref_scores value among the EM-selected references
    (selected_refs), or among all the contig's references if none of them was selected.
    """

    logging.info("Classifying contigs with thresholds KNOWN >= %.3f; NOVEL < %.3f.", known_thr, novel_thr)

    asm_lengths = read_assembly_lengths(assembly_fasta)

    # Keep PAF qlen where available and valid.
    for cid, paf_len in paf_contig_lengths.items():

        if paf_len > 0:
            asm_lengths[cid] = paf_len

    classification = {}

    counts = {"KNOWN": 0, "UNCERTAIN": 0, "NOVEL_CANDIDATE": 0}

    for cid, length in asm_lengths.items():

        if length < min_contig_len:
            continue

        scores = contig_scores.get(cid, {})

        if not scores:

            classification[cid] = ("NOVEL_CANDIDATE", 0.0, None)

            counts["NOVEL_CANDIDATE"] += 1

            continue

        best_score = max(scores.values())

        if ref_scores is None:
            best_ref = max(scores, key=scores.get)
        else:
            candidate_scores = ref_scores[cid]
            candidates = [ref for ref in candidate_scores if ref in selected_refs] or list(candidate_scores)
            best_ref = max(candidates, key=candidate_scores.get)

        if best_score >= known_thr:

            label = "KNOWN"

        elif best_score < novel_thr:

            label = "NOVEL_CANDIDATE"

        else:

            label = "UNCERTAIN"

        classification[cid] = (label, best_score, best_ref)

        counts[label] += 1

    logging.info("Contig classification: %d KNOWN | %d UNCERTAIN | " "%d NOVEL_CANDIDATE", counts["KNOWN"], counts["UNCERTAIN"], counts["NOVEL_CANDIDATE"])

    return classification, asm_lengths


# =============================================================================
# FASTA extraction
# =============================================================================


def extract_sequences(assembly_fasta, contig_ids, output_fasta):
    """
    Extract specified assembly contigs into a FASTA file.
    """

    wanted = set(contig_ids)

    written = set()

    current_id = None
    current_seq = []

    def flush(fh):

        if current_id is not None and current_id in wanted:

            fh.write(f">{current_id}\n")

            fh.write("".join(current_seq) + "\n")

            written.add(current_id)

    with open(assembly_fasta) as inp, open(output_fasta, "w") as out:

        for line in inp:

            line = line.rstrip()

            if line.startswith(">"):

                flush(out)

                current_id = norm_contig_id(line[1:].split()[0])

                current_seq = []

            else:

                current_seq.append(line)

        flush(out)

    missing = wanted - written

    if missing:

        logging.warning("%d requested extra contigs were not found in assembly.", len(missing))

    logging.info("Extracted %d extra contigs -> %s", len(written), output_fasta)

    return written


# =============================================================================
# Classification report
# =============================================================================


def write_classification_tsv(classification, contig_lengths, output_path):

    with open(output_path, "w") as fh:

        fh.write("contig_id\tlabel\tvcig_norm\tbest_ref\tcontig_length\n")

        for cid, (label, score, ref) in sorted(classification.items()):

            fh.write(f"{cid}\t" f"{label}\t" f"{score:.6f}\t" f"{ref or ''}\t" f"{contig_lengths.get(cid, 0)}\n")

    logging.info("Contig classification -> %s", output_path)


# =============================================================================
# Reduced database construction
# =============================================================================


def build_reduced_database(database, reduced_refs_path, reduced_db_path, threads):
    """
    Extract selected references from the full database using seqkit.
    """

    logging.info("Building reduced database.")

    cmd = ["seqkit", "grep", "-f", reduced_refs_path, "-j", str(threads), database]

    logging.info("Running: %s", " ".join(cmd))

    with open(reduced_db_path, "w") as output:

        subprocess.run(cmd, stdout=output, check=True)

    logging.info("Reduced database -> %s", reduced_db_path)


# =============================================================================
# Main module entry point
# =============================================================================


def run(args):
    """
    Run MADRe database reduction with genome-atomic reduced-reference
    selection.

    Returns
    -------
    dict
        classification
        contig_scores
        contig_lengths
        reduced_refs
        reduced_refs_path
        reduced_db_path
        extra_sequences
        classification_tsv
        mapping_class
        genome_expansion_added
    """

    os.makedirs(args.output_dir, exist_ok=True)

    logging.info("Starting DatabaseReduction.")

    logging.info("Assembly: %s", args.assembly)

    logging.info("Assembly mapping: %s", args.contig_paf)

    logging.info("Reference genome map: %s", args.reference_genome_map)

    # =========================================================================
    # 1. MADRe reduction
    # =========================================================================

    number_of_clusters = get_number_of_clusters(args.num_collapsed_strains)

    mapping_class_path = os.path.join(args.output_dir, "mapping_contig_class.txt")

    contig_vcig, paf_contig_lengths, references, contig_alignments = read_contig_paf(
        args.contig_paf, args.min_contig_len, no_mapq_filter=getattr(args, "no_mapq_filter", False)
    )

    U, NU, genomes_ids = calculate_mapping_class(contig_vcig, references, mapping_class_path=mapping_class_path)

    U, NU = clear_NU(U, NU, number_of_clusters, args.strictness, args.collapsed_strains_overhead)

    _, _, _, NU = EM(U, NU, list(genomes_ids.values()))

    reduced_refs = get_reduced_list(U, NU, genomes_ids, args.strain_species_json)

    # =========================================================================
    # 2. ND contig classification
    # =========================================================================

    novelty_score = getattr(args, "novelty_score", "union")
    logging.info(f"Novelty score: {novelty_score}")
    contig_scores = contig_novelty_scores(contig_vcig, paf_contig_lengths, contig_alignments, novelty_score)

    # union decides only the label; the best reference is taken from the EM-selected references (summed scores)
    ref_scores, selected_refs = None, None
    if novelty_score == "union":
        ref_scores = contig_novelty_scores(contig_vcig, paf_contig_lengths, method="summed")
        selected_refs = set(reduced_refs)

    classification, contig_lengths = classify_contigs(
        args.assembly, contig_scores, paf_contig_lengths, args.min_contig_len, args.known_thr, args.novel_thr, ref_scores, selected_refs
    )

    # =========================================================================
    # 3. Make sure confidently KNOWN references survive reduction
    # =========================================================================

    known_best_refs = {ref for label, _, ref in classification.values() if label == "KNOWN" and ref is not None}

    n_before = len(set(reduced_refs))

    reduced_refs = sorted(set(reduced_refs) | known_best_refs)

    n_added = len(reduced_refs) - n_before

    if n_added:

        logging.info("Added %d KNOWN best-reference(s) to reduced database.", n_added)

    # =========================================================================
    # 3.5. Genome-atomic expansion
    #
    # If any contig of a genome was selected above (by unique mapping,
    # species-aware NU selection, or as a confidently KNOWN best reference),
    # every other contig of that same genome is added too, so no selected
    # genome is left only partially represented in the reduced database.
    # =========================================================================

    ref_to_genome, genome_to_refs = load_reference_genome_map(args.reference_genome_map)

    n_before_expansion = len(reduced_refs)

    reduced_refs = expand_to_genome_level(reduced_refs, ref_to_genome, genome_to_refs)

    genome_expansion_added = len(reduced_refs) - n_before_expansion

    # =========================================================================
    # 4. Write reduced reference list
    # =========================================================================

    reduced_refs_path = os.path.join(args.output_dir, "reduced_db_refs.txt")

    write_reference_list(reduced_refs, reduced_refs_path)

    # =========================================================================
    # 5. Extract UNCERTAIN + NOVEL contigs
    # =========================================================================

    extra_ids = {cid for cid, (label, _, _) in classification.items() if label in {"UNCERTAIN", "NOVEL_CANDIDATE"}}

    extra_sequences_path = os.path.join(args.output_dir, "extra_sequences.fa")

    extract_sequences(args.assembly, extra_ids, extra_sequences_path)

    # =========================================================================
    # 6. Classification report
    # =========================================================================

    classification_tsv = os.path.join(args.output_dir, "contig_classification.tsv")

    write_classification_tsv(classification, contig_lengths, classification_tsv)

    # =========================================================================
    # 7. Optionally build actual reduced FASTA
    # =========================================================================

    reduced_db_path = None

    database = getattr(args, "database", None)

    if database is not None:

        reduced_db_path = os.path.join(args.output_dir, "reduced_db.fa")

        build_reduced_database(database, reduced_refs_path, reduced_db_path, getattr(args, "threads", 8))

    logging.info("DatabaseReduction complete.")

    return {
        "classification": classification,
        "contig_scores": contig_scores,
        "contig_lengths": contig_lengths,
        "reduced_refs": reduced_refs,
        "reduced_refs_path": reduced_refs_path,
        "reduced_db_path": reduced_db_path,
        "extra_sequences": extra_sequences_path,
        "classification_tsv": classification_tsv,
        "mapping_class": mapping_class_path,
        "genome_expansion_added": genome_expansion_added,
    }


# =============================================================================
# Standalone CLI
# =============================================================================


def main():

    parser = argparse.ArgumentParser(prog="DatabaseReduction.py", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    required = parser.add_argument_group("Required")

    required.add_argument("--assembly", required=True)

    required.add_argument("--contig_paf", required=True)

    required.add_argument("--strain_species_json", required=True)

    required.add_argument(
        "--reference_genome_map",
        required=True,
        help=(
            "reference_genome_map.tsv (reference, genome_id, "
            "genome_length), built by PrepareDatabase.py. Used to "
            "expand the reduced reference set to whole genomes."
        ),
    )

    required.add_argument("--num_collapsed_strains", required=True)

    required.add_argument("--output_dir", required=True)

    parser.add_argument("--database", default=None, help=("Full reference FASTA. If supplied, reduced_db.fa " "is built automatically."))

    parser.add_argument("--threads", type=int, default=8)

    parser.add_argument("--strictness", default="very-strict", choices=["less-strict", "strict", "very-strict"])

    parser.add_argument("--collapsed_strains_overhead", type=int, default=MAX_COLLAPSED_STRAINS_OVERHEAD)

    parser.add_argument("--known_thr", type=float, default=DEFAULT_KNOWN_THR)

    parser.add_argument("--novel_thr", type=float, default=DEFAULT_NOVEL_THR)

    parser.add_argument("--min_contig_len", type=int, default=MIN_CONTIG_LEN)

    parser.add_argument("--no_mapq_filter", action="store_true")

    parser.add_argument("--novelty-score", choices=["summed", "union"], default="union")

    run(parser.parse_args())


if __name__ == "__main__":
    main()
