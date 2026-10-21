# CalculateAbundances: abundance from StrainRefine's final read assignments (step 9 of MADRe).
#
# - Known species: genome copies per genome (reference_genome_map) = min(bases / genome length,
#   length-weighted median depth of the genome's sequences with reads); --genome_depth mean: bases / genome length.
#   Novel contigs placed in a species or higher taxon: bases / median genome length of that taxon.
# - Unplaced novel groups (NOVEL_COMPARE_ / UNRESOLVED_): bases / median genome length of the
#   group's candidate family, or of the whole database.
# - Host groups keep their reads and signal but are left out of the relative abundance.
# - Optional --min-abundance: species below it are left out of the taxonomic profile (all of them are kept in
#   strainrefine_taxonomic_profile.all_species.tsv).
# Outputs: strainrefine_reference_summary.tsv, strainrefine_species_or_group_summary.tsv and
# strainrefine_taxonomic_profile.tsv.
import logging
import math
import os
from statistics import median

import NovelPlacement
from NovelPlacement import ancestor, read_taxonomy
from StrainRefine import GROUP_PREFIXES, HOST_PREFIX, get_reference_id, load_reference_genome_map


def load_fasta_lengths(fasta_path):
    """
    Return {sequence_header: sequence_length} for a FASTA file.
    """
    lengths = {}

    if fasta_path is None:
        return lengths

    current = None
    current_len = 0

    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = current_len

                current = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())

    if current is not None:
        lengths[current] = current_len

    return lengths


def median_depth(sequence_bases, sequence_lengths):
    """
    Length-weighted median depth (bases / length) of a genome's sequences with reads: sort the sequences by
    depth and take the depth at which half of their summed length is reached. None if no length is known.
    """
    depths = sorted((bases / sequence_lengths[seq], sequence_lengths[seq]) for seq, bases in sequence_bases.items() if sequence_lengths.get(seq, 0) > 0)
    half, covered = sum(length for _, length in depths) / 2, 0
    for depth, length in depths:
        covered += length
        if covered >= half:
            return depth
    return None


def genome_copies(bases, length, sequence_bases=None, sequence_lengths=None):
    """
    Genome copies of one genome: bases / genome length (mean depth), or, with sequence_bases (sequence -> bases),
    min(mean depth, median_depth). The median never raises a genome; it lowers one whose extra reads sit on a small
    part of it, e.g. a high-copy plasmid.
    """
    mean = bases / length
    if not sequence_bases:
        return mean
    med = median_depth(sequence_bases, sequence_lengths or {})
    return mean if med is None else min(mean, med)


def calculate_known_genome_signals(
    group_refs, assigned_bases, genome_map, species_strain_info=None, diagnostics=None, taxonomy=None, genome_sequence_bases=None, sequence_lengths=None
):
    """
    Genome copies per known species (or GENUS_/TAXON_ group), summed over its genomes (genome_copies).
    Novel and host groups are handled in calculate_novel_group_signals.

    1. Every genome of the database's species (genome_map, also genomes without reads) -> its length;
       a GENUS_/TAXON_ group gets the genomes of all species under that taxon.
    2. References in genome_map: bases pooled per genome, signal = genome_copies (median depth when
       genome_sequence_bases, genome -> {sequence: bases}, is given; else bases / genome length).
       Other references (novel contigs placed in the group): bases / median genome length of the group.
    Returns signals (group -> genome copies) and denominators (summed length of the genomes with reads);
    diagnostics, if given, gets the parts of the signal for the species summary.
    """
    signals, denominators = {}, {}
    species_lengths = {}
    ref_groups = {ref: group for group, refs in group_refs.items() for ref in refs}
    for ref, (genome, length) in genome_map.items():
        group = species_strain_info.get(get_reference_id(ref)) if species_strain_info is not None else ref_groups.get(ref)
        if group is None or str(group).startswith(GROUP_PREFIXES):
            continue
        group = str(group)
        lengths = species_lengths.setdefault(group, {})
        if length <= 0 or (genome in lengths and lengths[genome] != length):
            raise ValueError(f"Invalid or inconsistent genome length for {genome}")
        lengths[genome] = length

    higher_groups = {g for g in group_refs if g.startswith(("GENUS_", "TAXON_"))}
    if higher_groups:
        if taxonomy is None:
            raise ValueError("Higher-rank abundance requires --taxonomy_nodes")
        parents, ranks = taxonomy
        for species, lengths in list(species_lengths.items()):
            if species.startswith(("GENUS_", "TAXON_")):
                continue
            lineage_ids = set()
            current = species
            while current not in lineage_ids:
                lineage_ids.add(current)
                if current not in parents:
                    break
                current = parents[current]
            for key in higher_groups:
                if key.split("_", 1)[1] in lineage_ids:
                    target = species_lengths.setdefault(key, {})
                    for genome, length in lengths.items():
                        if genome in target and target[genome] != length:
                            raise ValueError(f"Conflicting length for genome {genome}")
                        target[genome] = length

    for group, refs in group_refs.items():
        if group.startswith(GROUP_PREFIXES):
            continue
        genomes = {}
        unresolved_bases = 0
        unresolved_refs = 0
        for ref in refs:
            if ref not in genome_map or group.startswith(("GENUS_", "TAXON_")):
                unresolved_bases += assigned_bases.get(ref, 0)
                unresolved_refs += 1
                continue
            genome, length = genome_map[ref]
            bases, _ = genomes.get(genome, (0, length))
            genomes[genome] = (bases + assigned_bases.get(ref, 0), length)
        representative_length = 0
        estimated_signal = 0.0
        if unresolved_refs:
            lengths = species_lengths.get(group, {})
            if not lengths:
                raise ValueError(f"Species {group} has unresolved references but no mapped full genome length; extend --reference_genome_map.")
            representative_length = median(lengths.values())
            estimated_signal = unresolved_bases / representative_length
        seq_bases = genome_sequence_bases or {}
        mapped_signal = math.fsum(  # fsum: same result in any order
            genome_copies(bases, length, seq_bases.get(genome), sequence_lengths) for genome, (bases, length) in genomes.items()
        )
        signals[group] = mapped_signal + estimated_signal
        denominators[group] = sum(length for _, length in genomes.values())
        if diagnostics is not None:
            diagnostics[group] = {
                "mapped_genome_signal": mapped_signal,
                "unresolved_assigned_bases": unresolved_bases,
                "unresolved_reference_count": unresolved_refs,
                "representative_genome_length": representative_length,
                "estimated_genome_signal": estimated_signal,
                "estimated_signal_fraction": estimated_signal / signals[group] if signals[group] else 0.0,
            }
    return signals, denominators


def calculate_novel_group_signals(
    group_read_count, group_assigned_bases, group_final_refs, genome_map, taxonomy, contig_best_ref, ref_assigned_bases, diagnostics=None
):
    """
    Genome copies of unplaced novel groups (NOVEL_COMPARE_ / UNRESOLVED_): assigned bases / typical genome length.

    Typical length = median genome length of the group's candidate family (family of its contigs' best
    database hits, weighted by assigned bases), or the median of all genomes in genome_map if no family is known.
    The group's own contig length is not used: a few reads on a short contig would look like a whole genome.
    """
    parents, ranks = taxonomy
    novel_groups = [g for g in group_read_count if g.startswith(("NOVEL_COMPARE_", "UNRESOLVED_"))]
    if not novel_groups:
        return {}, {}

    family_cache = {}

    def family_of(taxid):
        if taxid not in family_cache:
            family_cache[taxid] = ancestor(taxid, parents, ranks, "family")
        return family_cache[taxid]

    group_family = {}
    for group in novel_groups:
        weight = {}
        for contig in group_final_refs.get(group, ()):
            best = contig_best_ref.get(contig, "")
            fam = family_of(get_reference_id(best)) if "|" in best else None
            if fam:
                weight[fam] = weight.get(fam, 0) + ref_assigned_bases.get(contig, 0)
        group_family[group] = max(sorted(weight), key=weight.get) if weight else None  # equal weight: smaller taxid

    wanted = {f for f in group_family.values() if f}
    genome_length = {}
    family_genomes = {}
    for ref, (genome, length) in genome_map.items():
        genome_length[genome] = length
        if wanted:
            fam = family_of(get_reference_id(ref))
            if fam in wanted:
                family_genomes.setdefault(fam, {})[genome] = length
    overall = median(genome_length.values()) if genome_length else 0

    signals, denominators = {}, {}
    for group in novel_groups:
        fam = group_family[group]
        lengths = family_genomes.get(fam)
        rep_len = median(lengths.values()) if lengths else overall
        bases = group_assigned_bases.get(group, 0)
        signals[group] = bases / rep_len if rep_len else 0.0
        denominators[group] = rep_len
        if diagnostics is not None:
            diagnostics[group] = {
                "mapped_genome_signal": 0.0,
                "unresolved_assigned_bases": bases,
                "unresolved_reference_count": len(group_final_refs.get(group, ())),
                "representative_genome_length": rep_len,
                "estimated_genome_signal": signals[group],
                "estimated_signal_fraction": 1.0,
            }
    logging.info(
        "Novel groups: %d unplaced groups counted as bases / typical genome length (%d with a candidate family, database median %.0f bp for the rest).",
        len(novel_groups),
        sum(1 for g in novel_groups if group_family[g] and family_genomes.get(group_family[g])),
        overall,
    )
    return signals, denominators


def calculate_abundance_data(
    final_assignments,
    species_strain_info,
    read_lengths,
    paf_reference_lengths,
    reference_fasta=None,
    reference_genome_map=None,
    taxonomy=None,
    contig_best_ref=None,
    read_sequences=None,
    genome_depth="median",
):
    """
    Reads, bases and genome-copy signal per reference and per species/group, and relative abundance.

    With reference_genome_map (always in MADRe):
        species, GENUS_, TAXON_                 calculate_known_genome_signals (bases / genome length)
        unplaced NOVEL_COMPARE_ / UNRESOLVED_   calculate_novel_group_signals (needs taxonomy and contig_best_ref)
    Groups not covered by these (host groups, novel groups without taxonomy, every group without a genome map)
    get the sequence signal: bases / summed length of the group's sequences.
    Relative abundance leaves host groups out of the total.
    genome_depth 'median' needs read_sequences (read -> the genome sequence it aligns to best, from StrainRefine);
    without it, or with 'mean', known genomes are bases / genome length.
    """
    reference_lengths = dict(paf_reference_lengths)
    if reference_fasta is not None:
        reference_lengths.update(load_fasta_lengths(reference_fasta))

    # reads and bases per final reference
    ref_read_count, ref_assigned_bases, missing_read_lengths = {}, {}, set()
    for read_id, final_ref in final_assignments:
        ref_read_count[final_ref] = ref_read_count.get(final_ref, 0) + 1
        read_len = read_lengths.get(read_id)
        if read_len is None:
            missing_read_lengths.add(read_id)
            continue
        ref_assigned_bases[final_ref] = ref_assigned_bases.get(final_ref, 0) + read_len
    if missing_read_lengths:
        logging.warning(
            "%d classified reads were not found in the PAF length table; their bases were not included in abundance calculations.", len(missing_read_lengths)
        )

    def group_of(ref):
        return species_strain_info.get(get_reference_id(ref))

    # reads, bases and references with reads per species/group
    group_read_count, group_assigned_bases, group_final_refs = {}, {}, {}
    for ref, count in ref_read_count.items():
        species_or_group = group_of(ref)
        if species_or_group is None:
            logging.warning("Reference %s (%s) not found in strain-species mapping; skipping it in abundance calculations.", ref, get_reference_id(ref))
            continue
        group = str(species_or_group)
        group_read_count[group] = group_read_count.get(group, 0) + count
        group_assigned_bases[group] = group_assigned_bases.get(group, 0) + ref_assigned_bases.get(ref, 0)
        group_final_refs.setdefault(group, set()).add(ref)

    # bases / own sequence length per reference (reported for references without a genome)
    ref_signal, missing_reference_lengths = {}, set()
    for ref, assigned_bases in ref_assigned_bases.items():
        ref_len = reference_lengths.get(ref)
        if ref_len is None or ref_len <= 0:
            missing_reference_lengths.add(ref)
            continue
        ref_signal[ref] = assigned_bases / ref_len
    if missing_reference_lengths:
        logging.warning("%d final references have no usable sequence length; their depth signal could not be calculated.", len(missing_reference_lengths))

    group_signal, group_denominator_length, genome_diagnostics = {}, {}, {}
    display_reference, display_reference_length, display_signal = {}, dict(reference_lengths), dict(ref_signal)
    best_database_hit = {}

    if reference_genome_map is not None:
        genome_map = load_reference_genome_map(reference_genome_map)
        genome_sequence_bases = None
        if genome_depth == "median" and read_sequences:
            # bases per sequence of each genome, each read on the sequence it aligns to best
            genome_sequence_bases = {}
            for read_id, final_ref in final_assignments:
                genome = genome_map.get(final_ref, (None, None))[0]
                if genome is None or read_id not in read_lengths:
                    continue
                seq = read_sequences.get(read_id, final_ref)
                per_seq = genome_sequence_bases.setdefault(genome, {})
                per_seq[seq] = per_seq.get(seq, 0) + read_lengths[read_id]
        elif genome_depth == "median":
            logging.warning("Median genome depth needs the per-read best sequences from StrainRefine; using bases / genome length.")
        known_signals, known_lengths = calculate_known_genome_signals(
            group_final_refs, ref_assigned_bases, genome_map, species_strain_info, genome_diagnostics, taxonomy, genome_sequence_bases, reference_lengths
        )
        if genome_sequence_bases:
            genome_length = dict(genome_map.values())
            lowered = sum(
                1
                for genome, per_seq in genome_sequence_bases.items()
                if genome_copies(sum(per_seq.values()), genome_length[genome], per_seq, reference_lengths)
                < 0.999 * sum(per_seq.values()) / genome_length[genome]
            )
            logging.info("Median genome depth: %d of %d genomes with reads lowered below bases / genome length.", lowered, len(genome_sequence_bases))
        group_signal.update(known_signals)
        group_denominator_length.update(known_lengths)
        if taxonomy is not None and contig_best_ref is not None:
            novel_signals, novel_lengths = calculate_novel_group_signals(
                group_read_count, group_assigned_bases, group_final_refs, genome_map, taxonomy, contig_best_ref, ref_assigned_bases, genome_diagnostics
            )
            group_signal.update(novel_signals)
            group_denominator_length.update(novel_lengths)
        # reference summary: a database reference is reported as its genome, with the whole genome length;
        # a contig gets the genome of its best database hit as a label only (best_database_hit, with its vcig)
        for ref in ref_read_count:
            genome_id, genome_length = genome_map.get(ref, (None, None))
            if genome_id is None:
                best = (contig_best_ref or {}).get(ref)
                inferred_genome_id, _ = genome_map.get(best, (None, None)) if best else (None, None)
                if inferred_genome_id is not None:
                    best_database_hit[ref] = inferred_genome_id
                continue
            display_reference[ref] = genome_id
            display_reference_length[ref] = genome_length
            display_signal[ref] = ref_assigned_bases.get(ref, 0) / genome_length if genome_length else 0.0
    else:
        if any(g.startswith(("GENUS_", "TAXON_")) for g in group_final_refs):
            raise ValueError("Genus abundance requires --reference_genome_map")
        logging.warning(
            "No reference genome map supplied: known-species abundance sums per-sequence depths and may be inflated for fragmented genomes. "
            "Use --reference_genome_map for genome-level normalization."
        )

    # sequence signal for the groups not covered above
    group_sequence_lengths = {}
    for ref, ref_len in reference_lengths.items():
        species_or_group = group_of(ref)
        if species_or_group is not None:
            group_sequence_lengths.setdefault(str(species_or_group), []).append(ref_len)
    for group in group_read_count:
        if group in group_signal:
            continue
        if group.startswith(GROUP_PREFIXES):
            total_length = sum(group_sequence_lengths.get(group, []))
            if total_length <= 0:
                logging.warning("Group %s has no usable sequence lengths; its abundance signal is set to 0.", group)
            group_signal[group] = group_assigned_bases.get(group, 0) / total_length if total_length > 0 else 0.0
            group_denominator_length[group] = total_length
        else:
            refs = group_final_refs.get(group, set())
            group_signal[group] = sum(ref_signal.get(ref, 0.0) for ref in refs)
            group_denominator_length[group] = sum(reference_lengths[ref] for ref in refs if ref in reference_lengths)

    # relative abundance of the microbial community: host groups keep their reads and signal but are left out of the total
    total_signal = math.fsum(signal for group, signal in group_signal.items() if not group.startswith(HOST_PREFIX))
    group_relative_abundance = {
        group: (100.0 * signal / total_signal if total_signal > 0 else 0.0) for group, signal in group_signal.items() if not group.startswith(HOST_PREFIX)
    }

    return {
        "genome_diagnostics": genome_diagnostics,
        "ref_read_count": ref_read_count,
        "ref_assigned_bases": ref_assigned_bases,
        "display_reference": display_reference,
        "display_reference_length": display_reference_length,
        "display_signal": display_signal,
        "best_database_hit": best_database_hit,
        "group_read_count": group_read_count,
        "group_assigned_bases": group_assigned_bases,
        "group_final_refs": group_final_refs,
        "group_denominator_length": group_denominator_length,
        "group_signal": group_signal,
        "group_relative_abundance": group_relative_abundance,
        "total_signal": total_signal,
    }


def write_reference_summary(abundance_data, ref_intervals, output_path):
    """
    strainrefine_reference_summary.tsv: one row per final reference with reads.
    A database reference is reported as its genome (genome_id, whole genome length, depth = bases / genome length);
    a contig keeps its own name and length, with the genome of its best database hit and the contig's vcig_norm
    (best_database_hit, best_database_hit_vcig: a label only, not an identification; low vcig = hardly similar).
    horizontal_coverage is the covered fraction of the reported sequence itself (not pooled per genome).
    """
    ref_count = abundance_data["ref_read_count"]
    display_reference = abundance_data["display_reference"]
    display_reference_length = abundance_data["display_reference_length"]
    display_signal = abundance_data["display_signal"]
    best_database_hit = abundance_data.get("best_database_hit", {})
    best_hit_vcig = abundance_data.get("best_database_hit_vcig", {})

    strain_notes = abundance_data.get("strain_notes", {})

    with open(output_path, "w") as f:

        f.write(
            "reference\t"
            "read_count\t"
            "assigned_bases\t"
            "reference_length\t"
            "depth_signal\t"
            "horizontal_coverage\t"
            "best_database_hit\t"
            "best_database_hit_vcig\t"
            "strain_note\n"
        )

        # most reads first; equal counts by reference name
        for ref, count in sorted(ref_count.items(), key=lambda x: (-x[1], x[0])):

            hcov = ref_intervals.get(ref, 0.0)

            assigned_bases = abundance_data["ref_assigned_bases"].get(ref, 0)

            reported_ref = display_reference.get(ref, ref)

            ref_len = display_reference_length.get(ref, 0)

            depth_signal = display_signal.get(ref, 0.0)

            hit_genome = best_database_hit.get(ref, "")

            hit_vcig = f"{best_hit_vcig[ref]:.6f}" if hit_genome and ref in best_hit_vcig else ""

            f.write(
                f"{reported_ref}\t"
                f"{count}\t"
                f"{assigned_bases}\t"
                f"{ref_len}\t"
                f"{depth_signal:.6f}\t"
                f"{hcov:.6f}\t"
                f"{hit_genome}\t"
                f"{hit_vcig}\t"
                f"{strain_notes.get(ref, '')}\n"
            )

    logging.info(f"Reference summary written to {output_path}")


def write_species_or_group_summary(abundance_data, output_path):
    """
    strainrefine_species_or_group_summary.tsv: one row per species or group (species, GENUS_, TAXON_, novel, host),
    highest abundance first. Host groups have relative_abundance NA.
    """
    group_counts = abundance_data["group_read_count"]
    group_assigned_bases = abundance_data["group_assigned_bases"]
    group_refs = abundance_data["group_final_refs"]
    group_denominator_length = abundance_data["group_denominator_length"]
    group_signal = abundance_data["group_signal"]
    group_relative_abundance = abundance_data["group_relative_abundance"]

    with open(output_path, "w") as f:

        f.write(
            "species_or_group\t"
            "read_count\t"
            "assigned_bases\t"
            "denominator_length\t"
            "genome_equivalent_signal\t"
            "relative_abundance\t"
            "num_references\t"
            "references\t"
            "mapped_genome_signal\tunresolved_assigned_bases\t"
            "unresolved_reference_count\trepresentative_genome_length\t"
            "estimated_genome_signal\testimated_signal_fraction\n"
        )

        # highest abundance first; equal abundance by name
        for species_or_group, count in sorted(group_counts.items(), key=lambda x: (-group_relative_abundance.get(x[0], 0.0), str(x[0]))):

            refs = sorted(group_refs.get(species_or_group, set()))

            rel_txt = "NA" if species_or_group.startswith(HOST_PREFIX) else f"{group_relative_abundance.get(species_or_group, 0.0):.6f}"
            diagnostic = abundance_data.get("genome_diagnostics", {}).get(species_or_group, {})
            diagnostic_values = "\t".join(
                str(diagnostic.get(key, "NA"))
                for key in (
                    "mapped_genome_signal",
                    "unresolved_assigned_bases",
                    "unresolved_reference_count",
                    "representative_genome_length",
                    "estimated_genome_signal",
                    "estimated_signal_fraction",
                )
            )
            f.write(
                f"{species_or_group}\t"
                f"{count}\t"
                f"{group_assigned_bases.get(species_or_group, 0)}\t"
                f"{group_denominator_length.get(species_or_group, 0)}\t"
                f"{group_signal.get(species_or_group, 0.0):.6f}\t"
                f"{rel_txt}\t"
                f"{len(refs)}\t"
                f"{','.join(refs)}\t{diagnostic_values}\n"
            )

    logging.info(f"Species/group summary written to {output_path}")


def load_ncbi_names(names_path):
    """
    Load scientific names from NCBI names.dmp.
    """
    names = {}

    with open(names_path, "r") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue

            taxid = parts[0]
            name_txt = parts[1]
            name_class = parts[3]

            if name_class == "scientific name":
                names[taxid] = name_txt

    return names


def strain_notes(updated_clusters, species_strain_info, contig_info, species_names, placed_species_contigs):
    """
    Strain note for the reference summary: an assembly contig in a known species that forms its own strain
    cluster (no database strain of that species in its cluster). It is below the KNOWN threshold to its
    best database hit, so it may be a new strain or a closely related new species. Information only:
    the species assignment is not changed.
    """
    notes = {}
    for cl in updated_clusters:
        if any("|" in ref for ref in cl):
            continue
        for ref in cl:
            group = str(species_strain_info.get(get_reference_id(ref), ""))
            if not group.isdigit() or ref not in contig_info:
                continue
            label, vcig, best = contig_info[ref]
            how = "placed by shared reads" if ref in placed_species_contigs else f"{label} contig"
            notes[ref] = (
                f"{how}; own strain cluster, no known strain of species {group} "
                f"({species_names.get(group, '?')}); "
                f"best hit vcig {vcig:.2f} (below the KNOWN threshold): new strain or "
                f"possibly a closely related new species"
            )
    if notes:
        logging.info("Strain notes: %d contigs flagged as divergent strains.", len(notes))
    return notes


def get_taxonomic_lineage(taxid, parents, ranks):
    """
    Return lineage entries from the supplied taxid upward.

    Only canonical ranks used for taxonomic profiling are retained.
    """
    wanted_ranks = {
        "acellular root",
        "cellular root",
        "other entries",
        "realm",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
        "superkingdom",
        "domain",
    }

    lineage = []
    seen = set()
    current = str(taxid)

    while current in parents and current not in seen:
        seen.add(current)

        rank = ranks.get(current, "no rank")

        if rank in wanted_ranks:
            normalized_rank = "domain" if rank == "superkingdom" else rank

            lineage.append((current, normalized_rank))

        parent = parents.get(current)

        if parent is None or parent == current:
            break

        current = parent

    return lineage


def write_taxonomic_profile(abundance_data, taxonomy, names, output_path, min_abundance=0.0):
    """
    strainrefine_taxonomic_profile.tsv: every species/GENUS_/TAXON_ signal added to all its NCBI ranks upward.
    Novel groups stay as rank 'novel_group' / 'unresolved'; host groups are left out.
    min_abundance > 0 (%): species rows below it are left out (no rescaling; higher ranks keep their full signal,
    as in the evaluation); the complete profile is also written to <output>.all_species.tsv.
    """

    parents, ranks = taxonomy

    group_counts = abundance_data["group_read_count"]
    group_signal = abundance_data["group_signal"]
    total_signal = abundance_data["total_signal"]

    taxon_counts = {}
    taxon_signals = {}
    taxon_ranks = {}

    for species_or_group, signal in group_signal.items():

        species_or_group_str = str(species_or_group)

        count = group_counts.get(species_or_group_str, 0)

        if species_or_group_str.startswith(HOST_PREFIX):
            continue

        if species_or_group_str.startswith(("NOVEL_COMPARE_", "UNRESOLVED_")):
            taxon_counts[species_or_group_str] = taxon_counts.get(species_or_group_str, 0) + count

            taxon_signals[species_or_group_str] = taxon_signals.get(species_or_group_str, 0.0) + signal

            taxon_ranks[species_or_group_str] = "unresolved" if species_or_group_str.startswith("UNRESOLVED_") else "novel_group"

            continue

        if species_or_group_str.startswith(("GENUS_", "TAXON_")):
            species_or_group_str = species_or_group_str.split("_", 1)[1]
            if ranks.get(species_or_group_str) not in {"genus", "family", "order", "class", "phylum", "domain", "superkingdom"}:
                raise ValueError(f"Invalid taxid in {species_or_group}: {species_or_group_str}")

        lineage = get_taxonomic_lineage(species_or_group_str, parents, ranks)

        if not lineage:
            if species_or_group_str not in parents:
                logging.warning("Taxid %s not found in nodes.dmp; writing it at species level only.", species_or_group_str)
            taxon_counts[species_or_group_str] = taxon_counts.get(species_or_group_str, 0) + count

            taxon_signals[species_or_group_str] = taxon_signals.get(species_or_group_str, 0.0) + signal

            taxon_ranks[species_or_group_str] = "species"

            continue

        for taxid, rank in lineage:
            taxon_counts[taxid] = taxon_counts.get(taxid, 0) + count

            taxon_signals[taxid] = taxon_signals.get(taxid, 0.0) + signal

            taxon_ranks[taxid] = rank

    rank_order = {
        "acellular root": 0,
        "cellular root": 0,
        "other entries": 0,
        "realm": 1,
        "domain": 1,
        "kingdom": 2,
        "phylum": 3,
        "class": 4,
        "order": 5,
        "family": 6,
        "genus": 7,
        "species": 8,
        "strain": 9,
        "novel_group": 10,
    }

    rows = []

    for taxid_or_group, signal in taxon_signals.items():

        rank = taxon_ranks.get(taxid_or_group, "unknown")

        if rank == "novel_group":
            name = taxid_or_group
        else:
            name = names.get(taxid_or_group, "")

        relative_abundance = 100.0 * signal / total_signal if total_signal > 0 else 0.0

        rows.append((rank_order.get(rank, 999), rank, -signal, str(taxid_or_group), name, taxon_counts.get(taxid_or_group, 0), signal, relative_abundance))

    rows.sort()

    header = "taxid_or_group\trank\tname\tread_count\tgenome_equivalent_signal\trelative_abundance\n"
    outputs = [(output_path, rows)]
    if min_abundance > 0:
        kept = [r for r in rows if not (r[1] == "species" and r[7] < min_abundance)]
        all_path = (output_path[:-4] if output_path.endswith(".tsv") else output_path) + ".all_species.tsv"
        outputs = [(output_path, kept), (all_path, rows)]
        dropped = [r for r in rows if r[1] == "species" and r[7] < min_abundance]
        logging.info(
            "Minimum abundance %.4g%%: %d of %d species left out of the profile (%.3f%% of the abundance); all species: %s",
            min_abundance,
            len(dropped),
            sum(1 for r in rows if r[1] == "species"),
            sum(r[7] for r in dropped),
            all_path,
        )
    for path, out_rows in outputs:
        with open(path, "w") as f:
            f.write(header)
            for _, rank, _, taxid_or_group, name, count, signal, relative_abundance in out_rows:
                f.write(f"{taxid_or_group}\t{rank}\t{name}\t{count}\t{signal:.6f}\t{relative_abundance:.6f}\n")

    logging.info(f"Taxonomic profile written to {output_path}")


def run(sr_result, args):
    """
    Abundance and the three abundance tables from StrainRefine.run()'s result.

    args: ref_summary_output, reference_genome_map, taxonomy_nodes, taxonomy_names, reference_fasta,
    contig_classification, and optionally genome_depth (median/mean, default median), min_abundance (%, default 0 = off),
    species_group_summary_output and taxonomic_profile_output (default: next to ref_summary_output).
    The returned abundance data also carries the taxonomy and names.
    """
    nodes = getattr(args, "taxonomy_nodes", None)
    names_path = getattr(args, "taxonomy_names", None)
    classification_path = getattr(args, "contig_classification", None)

    # each input file is read once
    taxonomy = read_taxonomy(nodes) if nodes else None
    names = load_ncbi_names(names_path) if names_path else {}
    contig_info = NovelPlacement.load_contig_classification(classification_path) if classification_path else {}
    contig_best_ref = {cid: best for cid, (_, _, best) in contig_info.items() if best} if classification_path else None

    abundance_data = calculate_abundance_data(
        sr_result["final_assignments"],
        sr_result["species_strain_info"],
        sr_result["read_lengths"],
        sr_result["reference_lengths"],
        reference_fasta=getattr(args, "reference_fasta", None),
        reference_genome_map=getattr(args, "reference_genome_map", None),
        taxonomy=taxonomy,
        contig_best_ref=contig_best_ref,
        read_sequences=sr_result.get("read_sequences"),
        genome_depth=getattr(args, "genome_depth", "median"),
    )
    abundance_data["strain_notes"] = strain_notes(
        sr_result["updated_clusters"], sr_result["species_strain_info"], contig_info, names, sr_result["placed_species_contigs"]
    )

    abundance_data["best_database_hit_vcig"] = {cid: vcig for cid, (_, vcig, best) in contig_info.items() if best}

    abundance_data["taxonomy"], abundance_data["names"] = taxonomy, names

    write_reference_summary(abundance_data, sr_result["ref_intervals"], args.ref_summary_output)
    write_outputs(abundance_data, args)

    return abundance_data


def write_outputs(abundance_data, args):
    """
    Species/group summary and taxonomic profile (the latter only with taxonomy and names), next to
    args.ref_summary_output unless species_group_summary_output / taxonomic_profile_output are given.
    """
    out_dir = os.path.dirname(args.ref_summary_output)
    species_group_summary_output = getattr(args, "species_group_summary_output", None) or os.path.join(out_dir, "strainrefine_species_or_group_summary.tsv")
    write_species_or_group_summary(abundance_data, species_group_summary_output)

    if abundance_data.get("taxonomy") is not None and getattr(args, "taxonomy_names", None) is not None:
        taxonomic_profile_output = getattr(args, "taxonomic_profile_output", None) or os.path.join(out_dir, "strainrefine_taxonomic_profile.tsv")
        write_taxonomic_profile(
            abundance_data, abundance_data["taxonomy"], abundance_data["names"], taxonomic_profile_output, getattr(args, "min_abundance", 0.0) or 0.0
        )
