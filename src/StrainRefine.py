# StrainRefine: read assignment and strain clustering (step 8 of MADRe).
# Abundance from its result is computed by CalculateAbundances (step 9).
#
# - Novel-contig placement (NovelPlacement.py, on by default): right after the
#   reads PAF is parsed, every NOVEL_COMPARE/UNRESOLVED contig is placed in the
#   taxonomy from the reads it shares with known references, and its group in
#   the species mapping becomes the species / GENUS_ / TAXON_ it is counted
#   at. A contig placed in a species takes part in that species' strain
#   clustering like any known reference (it joins a known strain's cluster or
#   forms its own = new strain). Outputs: novel_placement.tsv,
#   extended_strain_species.placed.json, novel_placement_report/ (dendrograms)
#   next to the read classification.
# - Genome pooling: before per-species reference clustering, the read-support
#   vectors of all references of one source genome (reference_genome_map.tsv)
#   are pooled into one vector, so a fragmented draft genome clusters as one
#   strain. See pool_reads_by_genome() and its call site in run().
import argparse
from collections import defaultdict
import logging
import json
from pathlib import Path
import random
import numpy as np
from sklearn.cluster import DBSCAN
from scipy.sparse import csr_matrix
import NovelPlacement
import os


# Host reference sequences (MADRe --host-reference) are named
# HOST|HOST_<name>|<original id>, so get_reference_id() gives HOST_<name>.
# Like novel groups, a host group is abundance-estimated from its own summed
# sequence length and is never looked up as an NCBI taxid.
HOST_PREFIX = "HOST_"
GROUP_PREFIXES = ("NOVEL_COMPARE_", "UNRESOLVED_", HOST_PREFIX)


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[logging.StreamHandler()])


def get_reference_id(reference):
    """
    Return the identifier used in the strain-species mapping.

    Standard MADRe references look like:
        |273121|NC_005090.1
    and use 273121 as the identifier.

    MADRe assembly-derived references look like:
        contig_68
    and use the full contig name as the identifier.
    """
    parts = reference.split("|")
    if len(parts) > 1 and parts[1]:
        return parts[1]
    return reference


def load_dict_from_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)


def species_split(U, NU, species_class, genomes, SS_info):
    species = set(species_class.values())

    genome_species = {ind: SS_info[get_reference_id(name)] for ind, name in genomes.items()}

    U_species = {}
    NU_species = {}

    for s in species:
        U_species[s] = {}
        NU_species[s] = {}

    for read_id, value_list in U.items():
        s = species_class[read_id]
        U_species[s][read_id] = value_list

    for read_id, value_list in NU.items():

        s = species_class[read_id]

        new_value_list = [[], [], [], 0]

        for i, ind in enumerate(value_list[0]):

            if genome_species[ind] == s:

                new_value_list[0].append(ind)

                new_value_list[1].append(value_list[1][i])

                new_value_list[2].append(value_list[2][i])

                new_value_list[3] = max(new_value_list[3], value_list[2][i])

        if not new_value_list[0]:

            logging.warning(f"Read {read_id} assigned to species {s} but has no mappings to it — skipping")

            continue

        NU_species[s][read_id] = new_value_list

    return U_species, NU_species


def species_identification_with_thresholds(U, NU, genomes, SS_info_json, min_read_count=5, min_mean_score=0.6, low_count_cap=30):

    logging.info("======Species identification======")

    logging.info(
        "Starting species identification with thresholds: "
        f"min_read_count={min_read_count}, "
        f"min_mean_score={min_mean_score}, "
        f"low_count_cap={low_count_cap}"
    )

    SS_info = load_dict_from_json(SS_info_json)

    all_mappings = {**U, **NU}

    logging.info(f"Processing {len(all_mappings)} reads ({len(U)} uniquely mapped, {len(NU)} ambiguously mapped)")

    def get_species_scores_for_read(read_id):

        species_scores = {}

        values = all_mappings[read_id]

        for i, genome_idx in enumerate(values[0]):

            genome = genomes[genome_idx]

            strain_taxid = get_reference_id(genome)

            species = SS_info[strain_taxid]

            score = values[2][i]

            if species not in species_scores:

                species_scores[species] = []

            species_scores[species].append((score, strain_taxid))

        return species_scores

    all_species_scores = {read_id: get_species_scores_for_read(read_id) for read_id in all_mappings}

    all_max_species_scores = {read_id: {sp: max(s[0] for s in scores) for sp, scores in ssd.items()} for read_id, ssd in all_species_scores.items()}

    all_sorted_species = {read_id: sorted(max_scores.items(), key=lambda x: x[1], reverse=True) for read_id, max_scores in all_max_species_scores.items()}

    del all_species_scores

    species_class = {}

    initial_ties = 0

    tie_rng = random.Random(0)  # fixed seed: tied reads are spread at random, the same way in every run

    for read_id in all_mappings.keys():

        max_scores = all_max_species_scores[read_id]

        best_species = None

        best_score = -1

        candidates = []

        for species, max_score in max_scores.items():

            if max_score > best_score:

                best_score = max_score

                best_species = species

                candidates = [species]

            elif max_score == best_score:

                candidates.append(species)

        if len(candidates) > 1:

            initial_ties += 1

            best_species = tie_rng.choice(candidates)

            logging.debug(f"Initial tie for read {read_id}: candidates={candidates}, chosen={best_species}")

        species_class[read_id] = best_species

    logging.info(f"Initial species assignment complete. Reads with tied best species score: {initial_ties}")

    stable = False

    max_iterations = 100

    iteration = 0

    # An iteration depends only on the assignments it starts from, so once the assignments repeat (after iteration i
    # they equal those after iteration i - period) every later iteration repeats too: (max_iterations - i) % period
    # more iterations give the assignments iteration max_iterations would end with. The last states are kept as
    # tuples in read order (the reads never change, only their species).
    max_period = 4

    recent_states = [tuple(species_class.values())]

    last_iteration = max_iterations

    while not stable and iteration < last_iteration:

        iteration += 1

        species_read_counts = {}

        species_score_sums = {}

        species_score_counts = {}

        for read_id, species in species_class.items():

            species_read_counts[species] = species_read_counts.get(species, 0) + 1

            max_scores = all_max_species_scores[read_id]

            if species in max_scores:

                species_score_sums[species] = species_score_sums.get(species, 0.0) + max_scores[species]

                species_score_counts[species] = species_score_counts.get(species, 0) + 1

        species_mean_scores = {s: (species_score_sums[s] / species_score_counts[s] if species_score_counts.get(s, 0) > 0 else 0.0) for s in species_read_counts}

        unreliable_species = set()

        for s in species_read_counts:

            if species_read_counts[s] < min_read_count or (species_mean_scores[s] < min_mean_score and species_read_counts[s] < low_count_cap):

                unreliable_species.add(s)

        logging.info(f"Iteration {iteration}: {len(species_read_counts)} active species, {len(unreliable_species)} unreliable species")

        if unreliable_species:

            logging.debug(
                f"Iteration {iteration} "
                "unreliable species: "
                + ", ".join([(f"{s}(count={species_read_counts[s]}, mean={species_mean_scores[s]:.4f})") for s in sorted(unreliable_species)])
            )

        reassigned_from_unreliable = 0

        reassigned_by_tie_break = 0

        for read_id in all_mappings.keys():

            current_species = species_class[read_id]

            if current_species in unreliable_species:

                for species, _ in all_sorted_species[read_id]:

                    if species not in unreliable_species:

                        if species_class[read_id] != species:

                            logging.debug(
                                f"Iteration {iteration}: "
                                f"read {read_id} reassigned "
                                f"from unreliable species "
                                f"{species_class[read_id]} "
                                f"to {species}"
                            )

                            reassigned_from_unreliable += 1

                        species_class[read_id] = species

                        break

            else:

                max_scores = all_max_species_scores[read_id]

                current_score = max_scores[current_species]

                best_alt_score = -1

                best_alt_species = None

                for species, sp_max in max_scores.items():

                    if species == current_species or species not in species_read_counts or species in unreliable_species:

                        continue

                    if sp_max > best_alt_score or (
                        sp_max == best_alt_score and species_read_counts.get(species, 0) > species_read_counts.get(best_alt_species, 0)
                    ):

                        best_alt_score = sp_max

                        best_alt_species = species

                if best_alt_species is not None:

                    if best_alt_score > current_score or (
                        best_alt_score == current_score and species_read_counts.get(best_alt_species, 0) > species_read_counts.get(current_species, 0)
                    ):

                        logging.debug(f"Iteration {iteration}: read {read_id} tie-resolved from {current_species} to {best_alt_species}")

                        species_class[read_id] = best_alt_species

                        reassigned_by_tie_break += 1

        changed_reads = reassigned_from_unreliable + reassigned_by_tie_break

        stable = changed_reads == 0

        logging.info(
            f"Iteration {iteration} summary: "
            f"{changed_reads} reads changed "
            "species assignment, "
            f"{reassigned_from_unreliable} "
            "reassigned from unreliable species, "
            f"{reassigned_by_tie_break} "
            "reassigned by tie-breaking, "
            f"stable={stable}"
        )

        if not stable and last_iteration == max_iterations:

            state = tuple(species_class.values())

            for period in range(2, max_period + 1):

                if len(recent_states) >= period and recent_states[-period] == state:

                    last_iteration = iteration + (max_iterations - iteration) % period

                    logging.info(
                        f"Iteration {iteration}: assignments repeat every {period} iterations; "
                        f"stopping after iteration {last_iteration}, which ends like iteration {max_iterations}"
                    )

                    break

            recent_states = (recent_states + [state])[-max_period:]

    final_species_counts = {}

    for species in species_class.values():

        final_species_counts[species] = final_species_counts.get(species, 0) + 1

    logging.info(f"Species identification finished after {iteration} iterations. Final number of assigned species: {len(final_species_counts)}")

    logging.debug("Final species counts: " + ", ".join([f"{s}:{c}" for s, c in sorted(final_species_counts.items(), key=lambda x: x[1], reverse=True)]))

    if not stable and last_iteration < max_iterations:

        logging.warning(f"Species identification did not converge (assignments repeat); result as after max_iterations={max_iterations}")

    elif not stable:

        logging.warning(f"Species identification reached max_iterations={max_iterations} before convergence")

    logging.info("==========================================")

    return species_class


def mapping_output(mapping_class_path, predictions_mapping):

    reduced = []

    if mapping_class_path is not None:

        f = open(mapping_class_path, "w")

        for contig, class_label in predictions_mapping.items():

            reduced.append(class_label)

            f.write(f"{contig} : {class_label}\n")

        f.close()

    logging.info("Mapping classification written in the file.")


def calculate_mapping_class(paf_path, species_strain_info, mapping_class_path=None, beta=2):

    logging.info("Mapping information extraction.")

    U = {}

    NU = {}

    read_genome_idx = {}

    nu_reads = set()

    genomes = {}

    genomes_names = {}

    genomes_id = 0

    ref_intervals = {}

    predictions_mapping = {}

    read_lengths = {}

    reference_lengths = {}

    with open(paf_path) as f:

        for line in f:

            parts = line.split()

            if len(parts) < 12:
                continue

            read_id = parts[0]

            read_lengths[read_id] = int(parts[1])

            reference_lengths[parts[5]] = int(parts[6])

            t_start = int(parts[7])

            t_end = int(parts[8])

            t_len = int(parts[6])

            length_t = t_end - t_start

            length_q = int(parts[3]) - int(parts[2])

            nm = int(parts[9])

            value_cig = nm / max(length_t, length_q)

            ref_prediction = parts[5]

            if get_reference_id(ref_prediction) not in species_strain_info:
                raise ValueError(f"Reference {ref_prediction} is not in the strain-species mapping")
            if ref_prediction not in ref_intervals:

                ref_intervals[ref_prediction] = ([], t_len)

            ref_intervals[ref_prediction][0].append((t_start, t_end))

            if ref_prediction not in genomes_names:

                genomes[genomes_id] = ref_prediction

                genomes_names[ref_prediction] = genomes_id

                genomes_id += 1

            ref_id = genomes_names[ref_prediction]

            if read_id not in read_genome_idx:

                U[read_id] = [[ref_id], [value_cig], [value_cig], value_cig]

                read_genome_idx[read_id] = {ref_id: 0}

                if mapping_class_path is not None:

                    predictions_mapping[read_id] = ref_prediction

            else:

                gidx = read_genome_idx[read_id]

                if ref_id in gidx:

                    idx = gidx[ref_id]

                    store = NU if read_id in nu_reads else U

                    if value_cig > store[read_id][1][idx]:

                        store[read_id][1][idx] = value_cig

                        store[read_id][2][idx] = value_cig

                        if value_cig > store[read_id][3]:

                            store[read_id][3] = value_cig

                else:

                    new_idx = len(gidx)

                    gidx[ref_id] = new_idx

                    if read_id not in nu_reads:

                        NU[read_id] = U.pop(read_id)

                        nu_reads.add(read_id)

                    NU[read_id][0].append(ref_id)

                    NU[read_id][1].append(value_cig)

                    NU[read_id][2].append(value_cig)

                    if value_cig > NU[read_id][3]:

                        NU[read_id][3] = value_cig

                    if mapping_class_path is not None and value_cig > NU[read_id][3]:

                        predictions_mapping[read_id] = ref_prediction

    if mapping_class_path is not None:

        mapping_output(mapping_class_path, predictions_mapping)

    for ref in ref_intervals:

        intervals, ref_len = ref_intervals[ref]

        if not intervals or ref_len == 0:

            ref_intervals[ref] = 0.0

            continue

        intervals.sort()

        cur_s, cur_e = intervals[0]

        covered = 0

        for s, e in intervals[1:]:

            if s <= cur_e:

                cur_e = max(cur_e, e)

            else:

                covered += cur_e - cur_s

                cur_s, cur_e = (s, e)

        covered += cur_e - cur_s

        ref_intervals[ref] = covered / ref_len

    # read and reference lengths for the abundance step, collected here so the PAF is read only once
    paf_info = {"read_lengths": read_lengths, "reference_lengths": reference_lengths}

    return (U, NU, genomes, ref_intervals, paf_info)


def pathoscope_redistribution(NU):

    for j in NU:

        scores = NU[j][2]

        xsum = sum(scores)

        if xsum == 0:

            NU[j][2] = [0.0] * len(scores)

        else:

            NU[j][2] = [s / xsum for s in scores]

    return NU


def find_medoid_and_avg_distance(cluster_indices, dist_matrix):

    if len(cluster_indices) == 1:

        return (cluster_indices[0], 0.0)

    sub_matrix = dist_matrix[np.ix_(cluster_indices, cluster_indices)]

    upper = sub_matrix[np.triu_indices_from(sub_matrix, k=1)]

    mean_distance = upper.mean()

    medoid_idx = cluster_indices[np.argmin(sub_matrix.mean(axis=1))]

    return (medoid_idx, mean_distance)


def initialize_clustering_output_dir(clustering_out):

    os.makedirs(clustering_out, exist_ok=True)

    return {
        "clusters_file": open(os.path.join(clustering_out, "clusters.txt"), "w"),
        "representatives_file": open(os.path.join(clustering_out, "representatives.txt"), "w"),
        "representatives_global": [],
        "avg_distances_global": [],
    }


def build_species_ref_dict(genomes, species_strain_info):

    species_ref_dict = {}

    references = set(genomes.values())

    for ref in references:

        taxid = get_reference_id(ref)

        species_taxid = species_strain_info[taxid]

        species_ref_dict.setdefault(species_taxid, []).append(ref)

    return species_ref_dict


def get_top_genomes(value_list, genomes):

    results = value_list[2]

    m = max(results)

    w = [i for i, x in enumerate(results) if x == m]

    p = [value_list[0][i] for i in w]

    return [genomes[i] for i in p]


def collect_species_read_data(all_mappings, genomes):

    ref_read_dict = {}

    reads_index_dict = {}

    classified = []

    reference_scores = {}

    ambigous_refs_count = {}

    ambigous_refs_reads = {}

    ambigous_species_count = 0

    ref_count = {}

    genomes_with_reads = set()

    relevant_genome_names = set(genomes[gidx] for vals in all_mappings.values() for gidx in vals[0])

    for name in relevant_genome_names:

        ref_read_dict[name] = set()  # positions of this species' reads that went to the reference (sparse read-support vector)

        ambigous_refs_count[name] = 0

        ambigous_refs_reads[name] = []

        ref_count[name] = 0

    for idx, read_id in enumerate(all_mappings.keys()):

        reads_index_dict[read_id] = idx

    for read_id, value_list in all_mappings.items():

        genome_list = get_top_genomes(value_list, genomes)

        if len(genome_list) == 1:

            genome = genome_list[0]

            classified.append((read_id, genome))

            reference_scores.setdefault(genome, [[], []])

            reference_scores[genome][0].append(read_id)

            reference_scores[genome][1].append(max(value_list[2]))

            ref_count[genome] += 1

            ambigous_refs_count[genome] += 1

            ref_read_dict[genome].add(reads_index_dict[read_id])

            genomes_with_reads.add(genome)

        else:

            for ref in genome_list:

                ambigous_refs_count[ref] += 1

                ambigous_refs_reads[ref].append(read_id)

                ref_read_dict[ref].add(reads_index_dict[read_id])

                genomes_with_reads.add(ref)

            ambigous_species_count += 1

    ref_read_dict = {ref: arr for ref in genomes_with_reads for arr in [ref_read_dict[ref]]}

    return {
        "ref_read_dict": ref_read_dict,
        "reads_index_dict": reads_index_dict,
        "classified": classified,
        "reference_scores": reference_scores,
        "ambigous_refs_count": ambigous_refs_count,
        "ambigous_refs_reads": ambigous_refs_reads,
        "ambigous_species_count": ambigous_species_count,
        "ref_count": ref_count,
    }


def resolve_ambiguous_reads(all_mappings, genomes, ambigous_refs_count, classified, reference_scores, ref_count):

    for read_id, value_list in all_mappings.items():

        read_scores = value_list[1]

        genome_list = get_top_genomes(value_list, genomes)

        if len(genome_list) <= 1:

            continue

        best_ref = genome_list[0]

        best_count = -1

        for ref in genome_list:

            c = ambigous_refs_count[ref]

            if c > best_count:

                best_ref = ref

                best_count = c

        classified.append((read_id, best_ref))

        ref_count[best_ref] += 1

        reference_scores.setdefault(best_ref, [[], []])

        reference_scores[best_ref][0].append(read_id)

        reference_scores[best_ref][1].append(max(read_scores))

    return (classified, reference_scores, ref_count)


def pool_reads_by_genome(ref_read_dict, reference_genome_map):
    """
    Before strain clustering, give all sequences of one genome the same read vector.

    ref_read_dict: reference sequence -> set of positions of this species' reads that went to this sequence.
    For every genome with two or more sequences (reference_genome_map), each sequence gets the union (OR) of
    the genome's read sets, so a draft genome's contigs, which see different reads, cluster as one strain.
    Sequences of different genomes are never combined; sequences missing from the map stay as they are.
    Only the clustering uses these vectors; read counts per sequence are not changed.
    """

    if not reference_genome_map:
        return ref_read_dict

    genome_members = defaultdict(list)

    for reference in ref_read_dict:

        genome_id = reference_genome_map.get(reference, (reference, None))[0]

        genome_members[genome_id].append(reference)

    pooled = dict(ref_read_dict)

    for references in genome_members.values():

        if len(references) < 2:
            continue

        union = set().union(*(ref_read_dict[reference] for reference in references))

        for reference in references:
            pooled[reference] = union

    return pooled


def jaccard_distance_matrix(read_sets):
    """
    Jaccard distance between every pair of read sets: 1 - shared / union.

    The sets are rows of a sparse 0/1 matrix M (references x reads), so the shared reads of all pairs are M @ M.T and
    memory grows with the number of read-reference assignments, not with references x reads. Same values as
    sklearn's pairwise_distances(dense 0/1 vectors, metric='jaccard') (same integer numerator and denominator).
    """
    indptr = np.cumsum([0] + [len(r) for r in read_sets])
    indices = np.fromiter((i for r in read_sets for i in r), dtype=np.int64, count=int(indptr[-1]))
    M = csr_matrix((np.ones(len(indices), dtype=np.int32), indices, indptr), shape=(len(read_sets), int(indices.max()) + 1 if len(indices) else 1))
    shared = (M @ M.T).toarray().astype(np.float64)
    size = np.diff(indptr).astype(np.float64)
    union = size[:, None] + size[None, :] - shared
    with np.errstate(invalid="ignore", divide="ignore"):
        dist = np.where(union > 0, (union - shared) / union, 0.0)
    np.fill_diagonal(dist, 0.0)
    return dist


def cluster_species_references(ref_read_dict, eps_value):

    if len(ref_read_dict) == 0:

        return {"filtered_ref_read_dict": {}, "clusters": [], "cluster_representatives": {}}

    ref_ids = list(ref_read_dict.keys())

    key_to_idx = {k: i for i, k in enumerate(ref_ids)}

    dist_matrix = jaccard_distance_matrix(list(ref_read_dict.values()))

    db = DBSCAN(metric='precomputed', eps=eps_value, min_samples=1)

    labels = db.fit_predict(dist_matrix)

    clusters = []

    cluster_representatives = {}

    for label in set(labels):

        cluster_indices = np.where(labels == label)[0]

        medoid_idx, _ = find_medoid_and_avg_distance(cluster_indices, dist_matrix)

        medoid = ref_ids[medoid_idx]

        cluster_refs = [ref_ids[i] for i in cluster_indices]

        ordered_cluster = [medoid] + [ref for ref in cluster_refs if ref != medoid]

        clusters.append(ordered_cluster)

        for ref in cluster_refs:

            cluster_representatives[ref] = medoid

    return {
        "filtered_ref_read_dict": ref_read_dict,
        "clusters": clusters,
        "cluster_representatives": cluster_representatives,
        "ref_ids": ref_ids,
        "key_to_idx": key_to_idx,
        "dist_matrix": dist_matrix,
    }


def summarize_cluster_support(clusters, reference_scores, high_score_threshold):

    ref_high_scores = {}

    ref_high_scores_global = {}

    low_score_read_count = 0

    assigned_reads = 0

    logging.info("======Cluster support summary======")

    for i, cluster in enumerate(clusters, start=1):

        logging.info(f"Cluster {i}:")

        for ref in cluster:

            if ref in reference_scores:

                scores = reference_scores[ref]

                high_scores = [s for s in scores[1] if (s >= high_score_threshold)]

                low_score_read_count += sum(1 for s in scores[1] if (s < high_score_threshold))

                ref_high_scores[ref] = len(high_scores)

                ref_high_scores_global[ref] = len(high_scores) / len(scores[1]) if len(scores[1]) > 0 else 0.0

                avg_score = sum(high_scores) / len(high_scores) if high_scores else 0.0

                assigned_reads += len(scores[0])

                logging.info(
                    "Reference: {}, "
                    "Assigned Reads: {}, "
                    "Average Score: {:.4f}, "
                    "High Scores: {}".format(ref, len(scores[0]), avg_score, len(high_scores))
                )

            else:

                ref_high_scores[ref] = 0

                ref_high_scores_global[ref] = 0.0

                logging.info("Reference: {}, Assigned Reads: 0, Average Score: 0.0000, High Scores: 0".format(ref))

    logging.info(f"Total low score reads globally: {low_score_read_count}")

    logging.info("==========================================")

    return {"ref_high_scores": ref_high_scores, "ref_high_scores_global": ref_high_scores_global, "assigned_reads": assigned_reads}


def choose_reference_reassignments(clusters, clusters_species, s_species_dist, ref_high_scores, ref_high_scores_global, cfg):

    changes = {}

    species_strong_cache = {}

    def is_weak(ref):
        high = ref_high_scores.get(ref, 0)
        fraction = ref_high_scores_global.get(ref, 0.0)
        return high < cfg["min_high_score_reads"] or (high < cfg["mid_high_score_reads"] and fraction < cfg["min_high_score_fraction"])

    # All weak references are known before any is moved, so a weak reference is never a target
    # and the result does not depend on the order in which clusters and references are visited.
    weak_refs = {ref for cluster in clusters for ref in cluster if is_weak(ref)}

    logging.info("======Reference reassignment summary======")

    for i, cluster in enumerate(clusters):

        species = clusters_species[i]

        (keys, key_to_idx, dist_matrix) = s_species_dist[species]

        if species not in species_strong_cache:

            strong_idx = np.array(
                [
                    idx
                    for idx, k in enumerate(keys)
                    if (ref_high_scores.get(k, 0) > cfg["min_high_score_reads"] or ref_high_scores_global.get(k, 0.0) > cfg["min_high_score_fraction"])
                ],
                dtype=int,
            )

            species_strong_cache[species] = (strong_idx, [keys[j] for j in strong_idx])

        (strong_idx, strong_keys) = species_strong_cache[species]

        for ref in cluster:

            if ref not in ref_high_scores:

                ref_high_scores[ref] = 0

            if ref not in ref_high_scores_global:

                ref_high_scores_global[ref] = 0.0

            weak_ref = ref_high_scores[ref] < cfg["min_high_score_reads"] or (
                ref_high_scores[ref] < cfg["mid_high_score_reads"] and ref_high_scores_global[ref] < cfg["min_high_score_fraction"]
            )

            if not weak_ref:

                continue

            ref_i = key_to_idx[ref]

            # closest strong reference first; equal distances by reference name
            for j in sorted(range(len(strong_idx)), key=lambda j: (dist_matrix[ref_i, strong_idx[j]], strong_keys[j])):

                cand = strong_keys[j]

                if cand == ref or cand in weak_refs:

                    continue

                changes[ref] = cand

                logging.info(
                    "Cluster {} - Reference {} "
                    "changed to {} based on "
                    "clustering with distance "
                    "{:.4f}".format(i + 1, ref, cand, dist_matrix[ref_i, strong_idx[j]])
                )

                break

    logging.info("==========================================")

    return changes


def get_ref_count(classified, changes):

    ref_count = {}

    for read_id, ref in classified:

        final_ref = changes[ref] if ref in changes else ref

        ref_count[final_ref] = ref_count.get(final_ref, 0) + 1

    return ref_count


def recompute_cluster_representatives_by_count(clusters, ref_count):

    changes = {}

    representatives_new = []

    updated_clusters = []

    for cluster in clusters:

        # most reads; equal counts by reference name
        max_ref = min(cluster, key=lambda r: (-ref_count.get(r, 0), r))

        if ref_count.get(max_ref, 0) > 0:

            for ref in cluster:

                changes[ref] = max_ref

            ordered = [max_ref] + [ref for ref in cluster if ref != max_ref]

        else:

            ordered = cluster

        representatives_new.append(ordered[0])

        updated_clusters.append(ordered)

    return (changes, representatives_new, updated_clusters)


def write_final_cluster_outputs(clustering_out, updated_clusters, representatives_new):

    with open(os.path.join(clustering_out, "clusters.txt"), "w") as f_clusters:

        # representative first, then the other members; clusters by representative (same file in every run)
        for cluster in sorted(updated_clusters, key=lambda c: c[0]):

            f_clusters.write(" ".join([cluster[0]] + sorted(cluster[1:])).strip() + "\n")

    with open(os.path.join(clustering_out, "representatives.txt"), "w") as f_rep:

        for rep in sorted(representatives_new):

            f_rep.write(rep + "\n")


def get_final_assignments(classified, changes, cluster_changes):
    """
    Final reference of every read: first the weak-reference reassignment (changes), then the
    strain-cluster representative (cluster_changes). Sorted by read id, so outputs are the same in every run.
    """
    final_assignments = []

    for read_id, ref in sorted(classified):

        ref = changes.get(ref, ref)

        final_assignments.append((read_id, cluster_changes.get(ref, ref)))

    return final_assignments


def best_genome_sequences(final_assignments, U, NU, genomes, reference_genome_map):
    """
    For every read finally assigned to a database genome: the sequence of that genome it aligns to best
    (highest vcig; equal vcig: the first in the PAF).

    Strain clustering gives all reads of a genome to one representative sequence, so per-sequence depth
    is lost in final_assignments; CalculateAbundances needs it for the median-depth genome copies.
    """
    if not reference_genome_map:
        return {}
    best = {}
    for read_id, ref in final_assignments:
        genome = reference_genome_map.get(ref, (None, None))[0]
        value = U.get(read_id) or NU.get(read_id)
        if genome is None or value is None:
            continue
        top, top_score = ref, -1.0
        for idx, score in zip(value[0], value[1]):
            seq = genomes[idx]
            if score > top_score and reference_genome_map.get(seq, (None, None))[0] == genome:
                top, top_score = seq, score
        best[read_id] = top
    return best


def write_final_assignments(final_assignments, output_path):

    with open(output_path, "w") as f:

        for read_id, final_ref in final_assignments:

            f.write(f"{read_id} : {final_ref}\n")


def load_reference_genome_map(path):
    """
    Read reference_genome_map.tsv (reference, genome_id, genome_length) -> {reference: (genome_id, genome_length)}.
    Stops on missing or inconsistent entries.
    """
    mapping = {}
    lengths = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        required = {"reference", "genome_id", "genome_length"}
        if not required.issubset(header):
            raise ValueError("Genome map requires reference, genome_id, genome_length columns")
        ref_col, genome_col, length_col = header.index("reference"), header.index("genome_id"), header.index("genome_length")
        n_cols = max(ref_col, genome_col, length_col) + 1
        for line_num, line in enumerate(handle, start=2):
            row = line.rstrip("\r\n").split("\t")
            if len(row) < n_cols:
                raise ValueError(f"Invalid genome map entry at line {line_num}")
            ref = row[ref_col].strip()
            genome = row[genome_col].strip()
            try:
                length = int(row[length_col])
            except ValueError:
                raise ValueError(f"Invalid genome length at line {line_num}") from None
            if not ref or not genome or length <= 0:
                raise ValueError(f"Invalid genome map entry at line {line_num}")
            if ref in mapping and mapping[ref] != (genome, length):
                raise ValueError(f"Conflicting genome mapping for {ref}")
            if genome in lengths and lengths[genome] != length:
                raise ValueError(f"Inconsistent full genome length for {genome}")
            mapping[ref] = (genome, length)
            lengths[genome] = length
    return mapping


def run(args):

    logging.info("Parameters:")

    logging.info(f"Strain-Species info JSON file path: {args.strain_species_info}")

    logging.info(f"Input PAF file path: {args.paf_path}")

    logging.info(f"Final classification labels file path: {args.read_class_output}")

    cfg = {
        "cluster_eps": args.cluster_eps,
        "high_score_threshold": args.high_score_threshold,
        "min_high_score_reads": args.min_high_score_reads,
        "mid_high_score_reads": args.mid_high_score_reads,
        "min_high_score_fraction": args.min_high_score_fraction,
        "species_min_read_count": args.species_min_read_count,
        "species_min_mean_score": args.species_min_mean_score,
        "species_low_count_cap": args.species_low_count_cap,
    }

    initialize_clustering_output_dir(args.clustering_out)

    species_strain_info = load_dict_from_json(args.strain_species_info)

    reference_genome_map_for_pooling = None

    if getattr(args, "reference_genome_map", None):

        reference_genome_map_for_pooling = load_reference_genome_map(args.reference_genome_map)

    (U, NU, genomes, ref_intervals, paf_info) = calculate_mapping_class(args.paf_path, species_strain_info, mapping_class_path=None, beta=0.5)

    placement = None
    if getattr(args, "novel_placement", True):
        if not (getattr(args, "taxonomy_nodes", None) and getattr(args, "contig_classification", None)):
            logging.warning("Novel-contig placement skipped: it needs --taxonomy_nodes and --contig_classification.")
        else:
            read_sets = defaultdict(set)
            for read_id, value in U.items():
                read_sets[genomes[value[0][0]]].add(read_id)
            for read_id, value in NU.items():
                for idx in value[0]:
                    read_sets[genomes[idx]].add(read_id)
            placement = NovelPlacement.Placement(
                read_sets,
                species_strain_info,
                NovelPlacement.load_contig_classification(args.contig_classification),
                args.taxonomy_nodes,
                workers=getattr(args, "threads", 1),
            )
            placement_changes = placement.compute()
            for contig, group in placement_changes.items():
                species_strain_info[contig] = group
            out_dir = Path(args.read_class_output).resolve().parent
            placed_json = out_dir / "extended_strain_species.placed.json"
            with open(placed_json, "w") as fh:
                json.dump(species_strain_info, fh)
            args.strain_species_info = str(placed_json)
            counted = {}
            for row in placement.rows:
                key = row.get("count_rank") or "not counted"
                counted[key] = counted.get(key, 0) + 1
            logging.info(
                "Novel-contig placement: %d contigs with a candidate family; counted at %s. Species mapping with placements: %s",
                len(placement.rows),
                counted,
                placed_json,
            )

    species_class = species_identification_with_thresholds(
        U,
        NU,
        genomes,
        args.strain_species_info,
        min_read_count=cfg["species_min_read_count"],
        min_mean_score=cfg["species_min_mean_score"],
        low_count_cap=cfg["species_low_count_cap"],
    )

    (U_species, NU_species) = species_split(U, NU, species_class, genomes, species_strain_info)

    species_ref_dict = build_species_ref_dict(genomes, species_strain_info)

    species = list(set(species_class.values()))

    s_species_dist = {}

    clusters = []

    clusters_species = []

    classified = []

    reference_scores = {}

    for s in species:

        logging.info("Species {} has {} references, eps value for clustering: {}".format(s, len(species_ref_dict[s]), cfg["cluster_eps"]))

        NU_result = pathoscope_redistribution(NU_species[s])

        all_mappings = {**U_species[s], **NU_result}

        species_data = collect_species_read_data(all_mappings, genomes)

        species_data["ref_read_dict"] = pool_reads_by_genome(species_data["ref_read_dict"], reference_genome_map_for_pooling)

        classified.extend(species_data["classified"])

        reference_scores.update(species_data["reference_scores"])

        (classified, reference_scores, _) = resolve_ambiguous_reads(
            all_mappings, genomes, species_data["ambigous_refs_count"], classified, reference_scores, species_data["ref_count"]
        )

        clustering_result = cluster_species_references(species_data["ref_read_dict"], cfg["cluster_eps"])

        if clustering_result["clusters"]:

            s_species_dist[s] = (clustering_result["ref_ids"], clustering_result["key_to_idx"], clustering_result["dist_matrix"])

        clusters.extend(clustering_result["clusters"])

        clusters_species.extend([s] * len(clustering_result["clusters"]))

    support_summary = summarize_cluster_support(clusters, reference_scores, cfg["high_score_threshold"])

    changes = choose_reference_reassignments(
        clusters, clusters_species, s_species_dist, support_summary["ref_high_scores"], support_summary["ref_high_scores_global"], cfg
    )

    ref_count = get_ref_count(classified, changes)

    (changes2, representatives_new, updated_clusters) = recompute_cluster_representatives_by_count(clusters, ref_count)

    active = [ref_count.get(rep, 0) > 0 for rep in representatives_new]

    updated_clusters = [c for c, keep in zip(updated_clusters, active) if keep]

    representatives_new = [r for r, keep in zip(representatives_new, active) if keep]

    write_final_cluster_outputs(args.clustering_out, updated_clusters, representatives_new)

    if placement is not None:
        placement.record_strain_outcome(updated_clusters)

    final_assignments = get_final_assignments(classified, changes, changes2)

    write_final_assignments(final_assignments, args.read_class_output)

    read_sequences = best_genome_sequences(final_assignments, U, NU, genomes, reference_genome_map_for_pooling)

    if placement is not None:
        final_reads = {}
        for _, ref in final_assignments:
            final_reads[ref] = final_reads.get(ref, 0) + 1
        out_dir = Path(args.read_class_output).resolve().parent
        placement.write(
            out_dir / "novel_placement.tsv",
            final_reads,
            report_dir=(out_dir / "novel_placement_report") if getattr(args, "taxonomy_names", None) and getattr(args, "placement_report", False) else None,
            names_path=getattr(args, "taxonomy_names", None),
        )
        logging.info("Novel-contig placement table: %s", out_dir / "novel_placement.tsv")

    # contigs NovelPlacement counted at species level (for the strain notes in CalculateAbundances)
    placed_species_contigs = {row["contig"] for row in placement.rows if row.get("count_rank") == "species"} if placement is not None else set()

    logging.info("Total assigned reads: {}".format(len(classified)))

    # everything CalculateAbundances needs
    return {
        "final_assignments": final_assignments,
        "read_sequences": read_sequences,
        "species_strain_info": species_strain_info,
        "read_lengths": paf_info["read_lengths"],
        "reference_lengths": paf_info["reference_lengths"],
        "ref_intervals": ref_intervals,
        "updated_clusters": updated_clusters,
        "placed_species_contigs": placed_species_contigs,
        "reference_genome_map": reference_genome_map_for_pooling,  # loaded once, reused by CalculateAbundances
        "taxonomy": (placement.parents, placement.ranks) if placement is not None else None,  # nodes.dmp, read by the placement
    }


def main():

    parser = argparse.ArgumentParser(description="MADRe.")

    parser.add_argument("--paf_path", type=str, required=True, help=("Path to the PAF file of assembly mapped to database."))

    parser.add_argument(
        "--reference_genome_map",
        default=None,
        help=(
            "TSV with reference, genome_id, genome_length columns. Exact FASTA "
            "reference IDs; genome_length is the full genome length in bp, "
            "including sequences removed by reduction. Pools known references "
            "within genomes before abundance normalization. Unmapped known "
            "references use pooled assigned bases divided by the species median "
            "full genome length (one length per mapped genome). Species without "
            "a mapped genome length raise an error. "
            "Unplaced NOVEL_COMPARE groups use assigned bases divided by the median "
            "genome length of their candidate family (or of the whole database). "
            "Without this option each reference's own sequence length is used."
        ),
    )

    parser.add_argument(
        "--no_novel_placement",
        dest="novel_placement",
        action="store_false",
        help=(
            "Do not place novel contigs in the taxonomy before species grouping "
            "(placement is on by default; it needs --taxonomy_nodes and "
            "--contig_classification)."
        ),
    )

    parser.add_argument(
        "--contig_classification",
        default=None,
        help=(
            "Optional path to DatabaseReduction's contig_classification.tsv "
            "(contig_id, label, vcig_norm, best_ref, contig_length columns). "
            "When set, a reference with no genome_id of its own in "
            "--reference_genome_map (e.g. a novel assembly contig chosen as "
            "a cluster representative) is labeled in the strain-level "
            "report's best_database_hit column (with its vcig_norm in best_database_hit_vcig) with the genome_id of "
            "its single best alignment hit, for readability only -- this "
            "never changes read counts, lengths, or abundance signal."
        ),
    )

    parser.add_argument(
        "--reference_fasta",
        type=str,
        default=None,
        help=(
            "FASTA used as the read-mapping target (normally the MADRe "
            "extended database). Recommended for abundance calculation so "
            "NOVEL_COMPARE group length includes contigs with zero PAF hits. "
            "If omitted, target lengths are taken from the PAF."
        ),
    )

    parser.add_argument("--strain_species_info", type=str, required=True, help=("JSON file with species/group information for every reference identifier."))

    parser.add_argument("--cluster_eps", type=float, default=0.8, help=("DBSCAN eps value for clustering references within species (default=0.8)."))

    parser.add_argument(
        "--high_score_threshold", type=float, default=0.6, help=("Score threshold used to count high-confidence assigned reads per reference (default=0.6).")
    )

    parser.add_argument(
        "--min_high_score_reads", type=int, default=5, help=("Minimum number of high-score reads for a reference to be considered supported (default=5).")
    )

    parser.add_argument("--mid_high_score_reads", type=int, default=10, help=("Intermediate support threshold for references (default=10)."))

    parser.add_argument(
        "--min_high_score_fraction", type=float, default=0.8, help=("Minimum fraction of high-score reads for moderately supported references (default=0.8).")
    )

    parser.add_argument("--species_min_read_count", type=int, default=5, help=("Minimum read count threshold for species reliability (default=5)."))

    parser.add_argument("--species_min_mean_score", type=float, default=0.6, help=("Minimum mean score threshold for species reliability (default=0.6)."))

    parser.add_argument(
        "--species_low_count_cap",
        type=int,
        default=30,
        help=("Species with mean score below threshold and read count below this cap are treated as unreliable (default=30)."),
    )

    parser.add_argument("--read_class_output", type=str, default="read_classification.out", help=("Path to output file with final classification labels."))

    parser.add_argument("--clustering_out", type=str, default="clustering_output", help=("Directory for clustering-related output files."))

    parser.add_argument(
        "--ref_summary_output", type=str, default="reference_summary.tsv", help=("TSV file with per-reference read counts and horizontal coverage.")
    )

    parser.add_argument(
        "--species_group_summary_output",
        type=str,
        default=None,
        help=(
            "TSV file with read counts aggregated by "
            "species or NOVEL_COMPARE group. "
            "If omitted, "
            "strainrefine_species_or_group_summary.tsv "
            "is written next to the reference summary."
        ),
    )

    parser.add_argument("--taxonomy_nodes", type=str, default=None, help=("Path to NCBI taxonomy nodes.dmp. Required for full taxonomic profile output."))

    parser.add_argument("--taxonomy_names", type=str, default=None, help=("Path to NCBI taxonomy names.dmp. Required for full taxonomic profile output."))

    parser.add_argument(
        "--taxonomic_profile_output",
        type=str,
        default=None,
        help=(
            "TSV file with abundance propagated across "
            "taxonomic ranks. If omitted, "
            "strainrefine_taxonomic_profile.tsv is written "
            "next to the reference summary."
        ),
    )

    parser.add_argument(
        "--min_abundance",
        type=float,
        default=0.0,
        help="Leave species below this relative abundance (%%) out of the taxonomic profile (complete profile: *.all_species.tsv). Default 0 (off).",
    )

    parser.add_argument(
        "--placement_report",
        action="store_true",
        help="Draw the per-contig figures in novel_placement_report/ (off by default; novel_placement.tsv is always written).",
    )
    # earlier option (the figures used to be drawn by default); accepted for old commands, no effect
    parser.add_argument("--no_placement_report", action="store_true", help=argparse.SUPPRESS)

    parser.add_argument("--threads", type=int, default=1, help="Processes for the novel-contig placement. Default: 1")

    parser.add_argument(
        "--genome_depth",
        choices=("median", "mean"),
        default="median",
        help=(
            "Genome copies of a known genome: 'median' (default) = min(bases / genome length, length-weighted median depth of its sequences "
            "with reads), so reads piled on a small high-copy sequence (plasmid, repeat) do not inflate the genome; 'mean' = bases / genome length."
        ),
    )

    args = parser.parse_args()

    result = run(args)

    import CalculateAbundances  # imported here: it imports helpers from this module

    CalculateAbundances.run(result, args)


if __name__ == "__main__":
    main()
