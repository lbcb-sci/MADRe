"""NovelPlacement: place novel assembly contigs in the taxonomy before
StrainRefine groups reads by species (MADRe / StrainRefine).

Everything is computed from StrainRefine's own read-to-reference data, and the
result is applied to the species mapping, so StrainRefine itself counts the
reads (a contig placed in a species also takes part in that species' strain
clustering).

For each NOVEL_COMPARE / UNRESOLVED contig with a candidate family (the family
of its best sequence hit, contig_classification.tsv), three agglomerative
clustering passes on shared-read sets (overlap-coefficient distance, average
linkage), finest first:
    pass 1  contig vs individual references of the family -> species
    pass 2  contig vs species (reads pooled per species)   -> genus
    pass 3  contig vs genera (reads pooled per genus)      -> family
then the counting rules:
    species kept only if the shared reads are >= species_min_contig_fraction
        of the contig's reads (else genus);
    genus kept only if the contig's best hit is in that genus with
        vcig >= genus_min_vcig (else family);
    sequence fallback: an unplaced contig, or (fallback_refines_family) a
        family-level one, gets its best hit's genus if vcig >= sequence_fallback_vcig.
The counted rank becomes the contig's group: species taxid, GENUS_<taxid>,
TAXON_<family taxid>, or unchanged (not counted).
"""

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform


def read_taxonomy(path):
    """NCBI nodes.dmp -> (parent of each taxid, rank of each taxid)."""
    parents, ranks = {}, {}
    with open(path) as handle:
        for line in handle:
            a = [v.strip() for v in line.split('|')]
            parents[a[0]], ranks[a[0]] = a[1], a[2]
    return parents, ranks


def ancestor(taxid, parents, ranks, rank='genus'):
    """The taxid itself or its first ancestor with the given rank; None if there is none."""
    seen = set()
    while taxid not in seen:
        seen.add(taxid)
        if ranks.get(taxid) == rank:
            return taxid
        if taxid not in parents:
            break
        taxid = parents[taxid]
    return None


NOVEL_PREFIXES = ('NOVEL_COMPARE_', 'UNRESOLVED_')


@dataclass
class PlacementOptions:
    cut_distance_pass1: float = 0.9
    cut_distance_pass2: float = 0.9
    cut_distance_pass3: float = 0.97
    min_population_size: int = 5
    min_shared_reads: int = 5
    min_contig_fraction: float = 0.0
    species_min_contig_fraction: float = 0.05
    genus_min_vcig: float = 0.1
    sequence_fallback_vcig: float = 0.2
    single_genus_family: bool = True
    fallback_refines_family: bool = True
    linkage_method: str = 'average'


def load_names(path):
    names = {}
    with open(path) as handle:
        for line in handle:
            a = [v.strip() for v in line.split('|')]
            if a[3] == 'scientific name':
                names[a[0]] = a[1]
    return names


def save_dendrogram(cid, pass_name, entities, z, cut_distance, names, out_path):
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    labels = []
    for taxid, reads in entities:
        if taxid == cid:
            labels.append(f'CONTIG {cid} ({len(reads)} reads)')
        else:
            labels.append(f'{names.get(taxid, taxid)} ({len(reads)} reads)')

    fig, ax = plt.subplots(figsize=(8, max(3, 0.4 * len(entities))))
    dendrogram(z, labels=labels, orientation='right', ax=ax)
    ax.axvline(cut_distance, color='red', linestyle='--', linewidth=1, label=f'cut distance = {cut_distance}')
    ax.set_xlabel('overlap-coefficient distance (0 = identical read-support, 1 = no overlap)')
    ax.set_title(f'{cid}: {pass_name}')
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# Genus label colours: first 7 slots of the validated categorical palette
# (dataviz reference palette; slot 8 red is reserved for the novel contig and
# its cluster). Genera beyond 7 fold into one grey "other genera" colour
# rather than cycling hues, which would make two genera look identical.
GENUS_PALETTE = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7']
OTHER_GENUS_COLOUR = '#7a7a7a'
CONTIG_COLOUR = '#e34948'
RANK_ORDER = {'species': 1, 'genus': 2, 'family': 3}
BRANCH_COLOUR = '#9a9a9a'


def overlap_distance(a, b):
    smaller = min(len(a), len(b))
    return 1.0 - len(a & b) / smaller if smaller else 1.0


def _draw_panel(ax, cid, contig_reads, shown, cut_distance, nm, genus_colour, title, colour_cluster=True):
    """One dendrogram panel. shown: [(taxon, reads, genus)] (contig excluded).
    Branches of the contig's cluster at cut_distance are drawn red."""
    entities = [(t, reads) for t, reads, _ in shown] + [(cid, contig_reads)]
    n = len(entities)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist[i, j] = dist[j, i] = overlap_distance(entities[i][1], entities[j][1])
    z = linkage(squareform(dist, checks=False), method='average')

    # Links of the subtree that forms the contig's flat cluster at the cut:
    # climb from the contig leaf while the merge height stays <= cut.
    children = {n + k: (int(z[k, 0]), int(z[k, 1])) for k in range(len(z))}
    parent = {c: node for node, cs in children.items() for c in cs}
    root, node = None, n - 1
    while node in parent and z[parent[node] - n, 2] <= cut_distance:
        node = parent[node]
        root = node
    red_links = set()
    stack = [root] if root is not None and colour_cluster else []
    while stack:
        k = stack.pop()
        if k >= n:
            red_links.add(k)
            stack.extend(children[k])

    labels = []
    for t, reads, g in shown:
        shared = len(contig_reads & reads)
        pct = 100 * shared / len(contig_reads) if contig_reads else 0
        tag = f' [{nm(g)}]' if g and g != t else ''
        labels.append(f'{nm(t)}{tag}  {len(reads)} reads, {shared} shared ({pct:.1f}% of contig)')
    labels.append(f'NOVEL {cid}  {len(contig_reads)} reads')

    dendrogram(z, labels=labels, orientation='right', ax=ax, link_color_func=lambda k: CONTIG_COLOUR if k in red_links else BRANCH_COLOUR)
    leaf_genus = {lab: g for lab, (_, _, g) in zip(labels, shown)}
    for tick in ax.get_yticklabels():
        text = tick.get_text()
        if text.startswith('NOVEL'):
            tick.set_color(CONTIG_COLOUR)
            tick.set_fontweight('bold')
        else:
            tick.set_color(genus_colour.get(leaf_genus.get(text), OTHER_GENUS_COLOUR))
    ax.axvline(cut_distance, color=CONTIG_COLOUR, linestyle='--', linewidth=1)
    ax.text(cut_distance + 0.006, 0.99, f'cut {cut_distance}', color=CONTIG_COLOUR, fontsize=8, transform=ax.get_xaxis_transform(), ha='left', va='top')
    ax.set_xlim(-0.02, 1.05)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_title(title, fontsize=9, loc='left', color='#333333', pad=8)
    in_cluster = root is not None
    return in_cluster


def write_contig_report(
    cid,
    fam,
    contig_reads,
    species_pools,
    decision,
    names,
    rank_of,
    cut_species,
    cut_genus,
    min_population_size,
    out_png,
    best_hit=None,
    best_vcig=0.0,
    pass_used='',
    min_shared_reads=5,
    min_contig_fraction=0.0,
    linkage_method='average',
    resolved=('', ''),
    final_reads=None,
    single_genus_family=False,
    ref_pools=None,
    species_of=None,
    cut_refs=0.9,
    rules=None,
    shared_fraction=None,
):
    """User-facing figure for one novel contig: which known species and genera
    (of its candidate family) its reads also map to, and how strongly.

    Panel 1: contig vs individual references of the family (pass 1, cut gives a species label; only the
    contig's cluster and its nearest references are drawn). Panel 2: contig vs known species (pooled per
    species, pass 2, cut gives a genus label). Panel 3: contig vs known genera (pooled per genus, pass 3,
    cut gives a family label). Red branches = the contig's cluster at that panel's cut. The panel that
    decided is marked, and the header names any rule that changed the rank afterwards.

    species_pools: {species_taxid: read_set}. Returns neighbour table rows."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    def nm(t):
        return names.get(t, t) if t else 'unknown'

    genus_pools = defaultdict(set)
    for sp, reads in species_pools.items():
        g = rank_of(sp, 'genus')
        if g:
            genus_pools[g] |= reads

    table = []
    for level, pools in (('species', species_pools), ('genus', genus_pools)):
        for t, reads in pools.items():
            shared = len(contig_reads & reads)
            g = rank_of(t, 'genus') if level == 'species' else t
            table.append(
                dict(
                    contig=cid,
                    level=level,
                    neighbour_taxid=t,
                    neighbour_name=nm(t),
                    neighbour_genus=nm(g),
                    distance=round(overlap_distance(contig_reads, reads), 4),
                    shared_reads=shared,
                    contig_reads=len(contig_reads),
                    neighbour_reads=len(reads),
                    fraction_of_contig_reads=round(shared / len(contig_reads), 4) if contig_reads else 0.0,
                    shown_in_figure=len(reads) >= min_population_size,
                )
            )
    table.sort(key=lambda r: (r['level'] != 'species', r['distance']))

    shown_sp = [(sp, reads, rank_of(sp, 'genus')) for sp, reads in species_pools.items() if len(reads) >= min_population_size]
    shown_g = [(g, reads, g) for g, reads in genus_pools.items() if len(reads) >= min_population_size]
    hidden = len(species_pools) - len(shown_sp)

    # Colour genera by closeness to the contig (most shared reads first), so
    # the fixed palette slots go to the genera that matter for this contig.
    genera = sorted(
        {g for _, _, g in shown_sp if g} | {g for g, _, _ in shown_g}, key=lambda g: (-len(contig_reads & genus_pools.get(g, set())), g)
    )  # equal shared reads: by taxid
    genus_colour = {g: (GENUS_PALETTE[i] if i < len(GENUS_PALETTE) else OTHER_GENUS_COLOUR) for i, g in enumerate(genera)}

    rank, label, source = decision
    n_final = final_reads if final_reads is not None else len(contig_reads)
    if rank:
        path = [(rank, label)] + [(r, rank_of(label, r)) for r in ('genus', 'family') if RANK_ORDER[r] > RANK_ORDER[rank]]
        path_txt = ' \u203a '.join(f'{r} {nm(t)}' for r, t in path if t)
        banner = f'{n_final} reads  \u2192  counted at {rank.upper()} {nm(label)}'
        if len([t for _, t in path if t]) > 1:
            banner += f'      ({path_txt})'
        rules = rules or {}
        hit_txt = nm(rank_of(best_hit, 'species') or best_hit) if best_hit else 'none'
        if source == 'clustering':
            why_line = 'decided by read clustering'
        elif source == 'sequence_refined':
            why_line = (f'clustering gave family {nm(resolved[1])}; refined to genus {nm(label)} by sequence: best hit {hit_txt} is in that family '
                        f'with vcig {best_vcig:.3f} >= {rules.get("sequence_fallback_vcig", 0.2)}')
        else:
            why_line = (f'no clustering pass gave a label; genus {nm(label)} from sequence alone: best hit {hit_txt} has vcig {best_vcig:.3f} '
                        f'>= {rules.get("sequence_fallback_vcig", 0.2)}')
    else:
        banner = f'{n_final} reads  \u2192  NOT COUNTED at any rank'
        why_line = (
            'no known species of this family has enough aligned reads to compare against'
            if not shown_sp
            else f'no clustering pass gave a label, and the sequence hit is too weak (vcig {best_vcig:.3f})'
        )
    r_rank, r_label = resolved
    if r_rank and (r_rank, r_label) != (rank, label) and source == 'clustering':
        rules = rules or {}
        if r_rank == 'species' and rank in ('genus', 'family'):
            frac = float(shared_fraction or 0)
            why_line += (f';  clustering gave species {nm(r_label)}, not kept: the shared reads are {100 * frac:.1f}% of the contig\'s reads '
                         f'(at least {100 * rules.get("species_min_contig_fraction", 0.05):.0f}% needed)')
        if rank == 'family':
            hit_genus = rank_of(best_hit, 'genus') if best_hit else None
            g_txt = nm(r_label) if r_rank == 'genus' else nm(rank_of(r_label, 'genus'))
            reason = (f'best hit {nm(rank_of(best_hit, "species") or best_hit)} has vcig {best_vcig:.3f} (at least {rules.get("genus_min_vcig", 0.1)} needed)'
                      if hit_genus and nm(hit_genus) == g_txt else f'best hit is in another genus ({nm(hit_genus)})')
            why_line += f';  genus {g_txt} not kept: a genus must be confirmed by the best sequence hit, and the {reason}'
    hit_sp = (rank_of(best_hit, 'species') or best_hit) if best_hit else None
    header = [
        f'{cid}: {len(contig_reads)} reads align to it (all alignments)   candidate family: {nm(fam)}',
        f'best sequence hit: {nm(hit_sp)} (vcig {best_vcig:.3f})',
    ]
    zero = sum(1 for reads in species_pools.values() if not reads)
    if hidden:
        few = hidden - zero
        parts = ([f'{zero} with no aligned reads'] if zero else []) + ([f'{few} with <{min_population_size} reads'] if few else [])
        header.append(f'not shown: {hidden} species of this family in the database ({", ".join(parts)})')

    def add_header(fig, top_in):
        """Banner + header lines in the top `top_in` inches of the figure."""
        h = fig.get_figheight()
        fig.text(
            0.01,
            1 - 0.12 / h,
            banner,
            fontsize=12,
            fontweight='bold',
            va='top',
            ha='left',
            color='white' if rank else '#222222',
            bbox=dict(boxstyle='round,pad=0.45', fc='#1f3b57' if rank else '#e3e3e3', ec='none'),
        )
        fig.text(0.01, 1 - 0.62 / h, why_line, fontsize=9, va='top', ha='left', color='#333333')
        fig.text(0.01, 1 - 0.92 / h, '\n'.join(header), fontsize=9, va='top', ha='left', color='#555555', linespacing=1.5)
        return 1 - top_in / h

    top_in = 1.05 + 0.24 * len(header)

    if not shown_sp:
        fig, ax = plt.subplots(figsize=(11, top_in + 0.8))
        ax.axis('off')
        ax.text(
            0.0,
            0.5,
            'No dendrogram: no known species of this family has enough aligned reads to ' 'compare against.',
            va='center',
            fontsize=9,
            color='#555555',
            transform=ax.transAxes,
        )
        fig.subplots_adjust(top=add_header(fig, top_in), bottom=0.02, left=0.01, right=0.99)
        fig.savefig(out_png, dpi=120, bbox_inches='tight')
        plt.close(fig)
        return table

    # Panel outcomes come from the same cluster_step the decision used, so
    # the picture cannot disagree with the result table.
    common = dict(
        min_population_size=min_population_size, min_shared_reads=min_shared_reads, min_contig_fraction=min_contig_fraction, linkage_method=linkage_method
    )
    why = {
        'too_few_entities': 'fewer than 2 known taxa to compare -> pass skipped',
        'alone_in_cluster': 'contig alone at the cut -> no label',
        'shared_reads_below_floor': f'cluster shares <{min_shared_reads} reads with contig -> no label',
        'contig_fraction_below_floor': 'cluster explains too little of the contig -> no label',
    }

    def outcome(cands, cut, level):
        members, _, status = cluster_step(cid, contig_reads, cands, cut_distance=cut, allow_single=single_genus_family and level == 'family', **common)
        if not members:
            return why.get(status, status), status != 'too_few_entities'
        if level == 'genus':
            gs = {rank_of(m, 'genus') for m in members}
            if len(gs) != 1 or None in gs:
                return f'cluster spans {len(gs)} genera -> no genus label', True
            return f'cluster = genus {nm(gs.pop())} -> genus label', True
        return f'cluster with {", ".join(nm(m) for m in members)} -> family {nm(fam)}', True

    used = str(pass_used)
    # pass 1: individual references; drawn = the contig's cluster plus its nearest other references
    shown_ref, t_ref, c_ref, n_ref_all = [], '', False, 0
    if ref_pools and species_of:
        def nm_ref(t):
            if '|' in t:
                sp = species_of(t.split('|')[1])
                return f'{nm(sp)} ({t.split("|")[-1]})'
            return nm(t)
        ref_cands = [(r, reads) for r, reads in ref_pools.items() if len(reads) >= min_population_size]
        n_ref_all = len(ref_cands)
        members, _, status = cluster_step(cid, contig_reads, ref_cands, cut_distance=cut_refs, **common)
        if members:
            sp = {species_of(m.split('|')[1]) for m in members}
            o_ref = f'cluster = species {nm(next(iter(sp)))} -> species label' if len(sp) == 1 else f'cluster spans {len(sp)} species -> no species label'
        else:
            o_ref = why.get(status, status)
        member_set = set(members or [])
        near = sorted((overlap_distance(contig_reads, reads), r) for r, reads in ref_cands if r not in member_set)
        keep = sorted(member_set, key=lambda r: overlap_distance(contig_reads, ref_pools[r]))[:20] + [r for _, r in near[:12]]
        shown_ref = [(r, ref_pools[r], rank_of(species_of(r.split('|')[1]), 'genus')) for r in keep]
        c_ref = bool(members)
        t_ref = (f'Reference view (pass 1, cut {cut_refs}): {o_ref}' + ('   [DECIDED HERE]' if used == '1' else '')
                 + (f'   (showing {len(shown_ref)} of {n_ref_all} references: its cluster and the nearest others)' if n_ref_all > len(shown_ref) else ''))
    o_sp, c_sp = outcome(list(species_pools.items()), cut_species, 'genus')
    o_g, c_g = outcome(list(genus_pools.items()), cut_genus, 'family')
    t_sp = f'Species view (pass 2, cut {cut_species}): {o_sp}' + ('   [DECIDED HERE]' if used == '2' else '') + ('   (not needed: decided in pass 1)' if used == '1' else '')
    t_g = (f'Genus view (pass 3, cut {cut_genus}): {o_g}' + ('   [DECIDED HERE]' if used == '3' else '')
           + ('   (not needed: decided in an earlier pass)' if used in ('1', '2') else ''))
    n_ref, n_sp, n_g = len(shown_ref) + 1, len(shown_sp) + 1, len(shown_g) + 1
    panels = ([('ref', n_ref)] if shown_ref else []) + [('sp', n_sp), ('g', n_g)]
    fig, axes = plt.subplots(len(panels), 1, figsize=(14, 0.38 * sum(n for _, n in panels) + 1.3 * len(panels) + top_in),
                             gridspec_kw=dict(height_ratios=[n + 1 for _, n in panels]))
    ax = dict(zip([k for k, _ in panels], axes))
    if shown_ref:
        _draw_panel(ax['ref'], cid, contig_reads, shown_ref, cut_refs, nm_ref, genus_colour, t_ref, c_ref)
    _draw_panel(ax['sp'], cid, contig_reads, shown_sp, cut_species, nm, genus_colour, t_sp, c_sp)
    if len(shown_g) >= 1:
        _draw_panel(ax['g'], cid, contig_reads, shown_g, cut_genus, nm, genus_colour, t_g, c_g)
    else:
        ax['g'].axis('off')
    ax['g'].set_xlabel('distance: 1 - shared reads / reads of smaller set   (0 = same reads, 1 = none shared)')

    handles = [Patch(color=genus_colour[g], label=nm(g)) for g in genera[: len(GENUS_PALETTE)]]
    if len(genera) > len(GENUS_PALETTE):
        handles.append(Patch(color=OTHER_GENUS_COLOUR, label=f'other genera ({len(genera) - len(GENUS_PALETTE)})'))
    handles.append(Patch(color=CONTIG_COLOUR, label="novel contig / its cluster at the cut"))
    fig.legend(handles=handles, loc='lower center', ncol=min(len(handles), 5), fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.04, 1, add_header(fig, top_in)))
    fig.savefig(out_png, dpi=120, bbox_inches='tight')
    plt.close(fig)
    return table


def candidate_distances(candidates, min_population_size):
    """Distances among the candidates cluster_step keeps (read set of at least min_population_size reads), in its order:
    1 - shared reads / reads of the smaller set (1 if either is empty), from a sparse read-membership matrix product.
    Same numbers as cluster_step's pairwise loop (exact integer counts, the same float division); computed once per
    family and pass instead of once per contig."""
    entities = [(label, reads) for label, reads in candidates if len(reads) >= min_population_size]
    labels = [label for label, _ in entities]
    n = len(entities)
    if n == 0:
        return labels, np.zeros((0, 0))
    index = {}
    rows, cols = [], []
    for i, (_, reads) in enumerate(entities):
        for r in reads:
            rows.append(i)
            cols.append(index.setdefault(r, len(index)))
    member = csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)), shape=(n, len(index)))
    shared = (member @ member.T).toarray().astype(np.float64)
    sizes = np.array([len(reads) for _, reads in entities], dtype=np.float64)
    smaller = np.minimum.outer(sizes, sizes)
    with np.errstate(divide='ignore', invalid='ignore'):
        dist = np.where(smaller > 0, 1.0 - shared / smaller, 1.0)
    np.fill_diagonal(dist, 0.0)
    return labels, dist


def cluster_step(
    cid,
    contig_raw_reads,
    candidates,
    min_population_size,
    min_shared_reads,
    min_contig_fraction,
    linkage_method,
    cut_distance,
    plot_info=None,
    allow_single=False,
    known=None,
):
    """candidates: [(label, read_set), ...], not including the contig itself.

    Returns (members, shared_count, status). members is a list of candidate
    labels sharing the contig's cluster, or None if no usable match.
    """
    entities = [(label, reads) for label, reads in candidates if len(reads) >= min_population_size]
    entities.append((cid, contig_raw_reads))
    n = len(entities)
    if n == 2 and allow_single:
        # Only one known taxon to compare with (e.g. a family with a single
        # known genus): no tree to cut, so test the contig against it directly.
        label, reads = entities[0]
        d = overlap_distance(contig_raw_reads, reads)
        shared = len(contig_raw_reads & reads)
        if d > cut_distance:
            return None, shared, 'alone_in_cluster'
        if shared < min_shared_reads:
            return None, shared, 'shared_reads_below_floor'
        if shared / len(contig_raw_reads) < min_contig_fraction:
            return None, shared, 'contig_fraction_below_floor'
        return [label], shared, 'ok'
    if n < 3:
        return None, 0, 'too_few_entities'

    if known is not None and known[0] == [label for label, _ in entities[:-1]]:
        # distances among the candidates computed once (candidate_distances); only the contig's own row here
        dist = np.zeros((n, n))
        dist[: n - 1, : n - 1] = known[1]
        for i in range(n - 1):
            a = entities[i][1]
            smaller = min(len(a), len(contig_raw_reads))
            dist[i, n - 1] = dist[n - 1, i] = 1.0 - (len(a & contig_raw_reads) / smaller) if smaller else 1.0
    else:
        dist = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                a, b = entities[i][1], entities[j][1]
                smaller = min(len(a), len(b))
                d_ij = 1.0 - (len(a & b) / smaller) if smaller else 1.0
                dist[i, j] = dist[j, i] = d_ij

    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method=linkage_method)
    labels = fcluster(z, t=cut_distance, criterion='distance')

    if plot_info is not None:
        pass_name, names, plot_path = plot_info
        save_dendrogram(cid, pass_name, entities, z, cut_distance, names, plot_path)

    contig_idx = n - 1
    contig_cluster = labels[contig_idx]
    members = [entities[i][0] for i in range(n) if labels[i] == contig_cluster and i != contig_idx]
    if not members:
        return None, 0, 'alone_in_cluster'

    member_reads = set()
    for i in range(n):
        if labels[i] == contig_cluster and i != contig_idx:
            member_reads |= entities[i][1]
    shared = len(contig_raw_reads & member_reads)
    if shared < min_shared_reads:
        return None, shared, 'shared_reads_below_floor'
    # The distance is normalized by the smaller set, which is often a tiny
    # known species (a handful of stray reads): 11 shared reads out of its 27
    # looks like 41% overlap while explaining 1% of the contig. Require the
    # cluster to explain a real share of the contig's own reads.
    if shared / len(contig_raw_reads) < min_contig_fraction:
        return None, shared, 'contig_fraction_below_floor'

    return members, shared, 'ok'


def load_contig_classification(path):
    """contig_id -> (label, vcig_norm, best_ref)."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            ref = (row.get('best_ref') or '').strip()
            out[row['contig_id']] = (row.get('label', ''), float(row.get('vcig_norm') or 0), ref if ref and ref.lower() != 'none' else '')
    return out


class Placement:
    """
    Novel-contig placement for one run, used by StrainRefine in this order:
    compute() -> {contig: new group}, which StrainRefine writes into its species mapping;
    after strain clustering, record_strain_outcome(); finally write() for novel_placement.tsv and the report.
    """

    def __init__(self, read_sets, species_info, classification, nodes_path, opts=None):
        self.opts = opts or PlacementOptions()
        self.read_sets = read_sets
        self.species_info = species_info
        self.classification = classification
        self.parents, self.ranks = read_taxonomy(nodes_path)
        self._rank_cache = {}
        self.rows = []
        self.contig_family = {}
        self.family_species_refs = defaultdict(lambda: defaultdict(set))
        self._species_pool_cache = {}
        self._known_cache = {}  # (pass, family) -> distances among that pass's candidates (candidate_distances)

    def rank_of(self, taxid, rank):
        key = (taxid, rank)
        if key not in self._rank_cache:
            self._rank_cache[key] = ancestor(taxid, self.parents, self.ranks, rank)
        return self._rank_cache[key]

    def species_of(self, taxid):
        return self.rank_of(taxid, 'species') or taxid

    def _known(self, key, candidates):
        """Distances among one pass's candidates of a family, computed on first use (the same for every contig of the family)."""
        if key not in self._known_cache:
            self._known_cache[key] = candidate_distances(candidates, self.opts.min_population_size)
        return self._known_cache[key]

    def species_pools(self, fam):
        """species taxid -> reads of all its references, for the species of family fam
        (built once per family; used by pass 2, pass 3 and the report)."""
        if fam not in self._species_pool_cache:
            pools = {}
            for sp_taxid, refs in self.family_species_refs.get(fam, {}).items():
                pooled = set()
                for r in refs:
                    pooled |= self.read_sets.get(r, set())
                pools[sp_taxid] = pooled
            self._species_pool_cache[fam] = pools
        return self._species_pool_cache[fam]

    def _count_rank(self, cid, rank, label, shared_fraction):
        o = self.opts
        _, vcig, ref = self.classification.get(cid, ('', 0.0, ''))
        if rank == 'species' and shared_fraction < o.species_min_contig_fraction:
            rank, label = 'genus', self.rank_of(label, 'genus')
        if rank == 'genus':
            hit_genus = self.rank_of(ref.split('|')[1], 'genus') if '|' in ref else None
            if not (label and hit_genus == label and vcig >= o.genus_min_vcig):
                rank, label = 'family', self.rank_of(label, 'family') if label else None
        return (rank, label) if label else ('', '')

    def compute(self):
        """Run the passes for every novel contig. Returns {contig: new_group}."""
        o, rank_of, species_of = self.opts, self.rank_of, self.species_of
        novel = {k: str(v) for k, v in self.species_info.items() if isinstance(v, str) and v.startswith(NOVEL_PREFIXES)}
        for cid in novel:
            ref = self.classification.get(cid, ('', 0.0, ''))[2]
            if '|' in ref:
                fam = rank_of(ref.split('|')[1], 'family')
                if fam:
                    self.contig_family[cid] = fam
        family_refs = defaultdict(list)
        for ref in self.read_sets:
            if '|' not in ref:
                continue
            taxid = ref.split('|')[1]
            fam = rank_of(taxid, 'family')
            if fam:
                family_refs[fam].append(ref)
                self.family_species_refs[fam][species_of(taxid)].add(ref)
        self.family_refs = family_refs
        common = dict(
            min_population_size=o.min_population_size,
            min_shared_reads=o.min_shared_reads,
            min_contig_fraction=o.min_contig_fraction,
            linkage_method=o.linkage_method,
        )
        for cid, fam in self.contig_family.items():
            contig_reads = self.read_sets.get(cid, set())
            _, vcig, best = self.classification.get(cid, ('', 0.0, ''))
            base = dict(
                contig=cid, original_group=novel[cid], candidate_family=fam, contig_reads=len(contig_reads), best_hit=best, best_hit_vcig=round(vcig, 4)
            )
            if len(contig_reads) < o.min_population_size:
                self.rows.append(dict(base, decision='contig_population_too_small'))
                continue
            result, status_log = None, []
            cands = [(r, self.read_sets.get(r, set())) for r in family_refs.get(fam, [])]
            members, shared, status = cluster_step(
                cid, contig_reads, cands, cut_distance=o.cut_distance_pass1, known=self._known(('references', fam), cands), **common
            )
            status_log.append(f'pass1={status}')
            if members:
                sp = {species_of(m.split('|')[1]) for m in members}
                if len(sp) == 1:
                    result = ('species', sp.pop(), members, shared, 1)
            if result is None:
                cands = list(self.species_pools(fam).items())
                members, shared, status = cluster_step(
                    cid, contig_reads, cands, cut_distance=o.cut_distance_pass2, known=self._known(('species', fam), cands), **common
                )
                status_log.append(f'pass2={status}')
                if members:
                    genera = {rank_of(m, 'genus') for m in members}
                    if len(genera) == 1 and None not in genera:
                        result = ('genus', genera.pop(), members, shared, 2)
            if result is None:
                pools = defaultdict(set)
                for sp_taxid, pooled in self.species_pools(fam).items():
                    g = rank_of(sp_taxid, 'genus')
                    if g:
                        pools[g] |= pooled
                members, shared, status = cluster_step(
                    cid,
                    contig_reads,
                    list(pools.items()),
                    cut_distance=o.cut_distance_pass3,
                    allow_single=o.single_genus_family,
                    known=self._known(('genera', fam), list(pools.items())),
                    **common,
                )
                status_log.append(f'pass3={status}')
                if members:
                    result = ('family', fam, members, shared, 3)
            if result is None:
                row = dict(base, decision='no_match_at_any_level (' + ', '.join(status_log) + ')')
            else:
                rank, label, members, shared, pass_num = result
                frac = shared / len(contig_reads)
                c_rank, c_label = self._count_rank(cid, rank, label, frac)
                row = dict(
                    base,
                    resolved_rank=rank,
                    resolved_label=label,
                    pass_used=pass_num,
                    cluster_members=','.join(str(m) for m in members),
                    contig_shared_fraction=round(frac, 4),
                    count_rank=c_rank,
                    count_label=c_label,
                    count_source='clustering' if c_rank else '',
                    decision='clustered',
                )
            genus = rank_of(best.split('|')[1], 'genus') if '|' in best else None
            if row.get('count_rank'):
                # only when clustering itself stopped at family; a genus or species from clustering that
                # _count_rank demoted (best hit in another genus) is not overruled by the best hit
                if (
                    o.fallback_refines_family
                    and row['count_rank'] == 'family'
                    and row['resolved_rank'] == 'family'
                    and genus
                    and rank_of(genus, 'family') == row['count_label']
                    and vcig >= o.sequence_fallback_vcig
                ):
                    row.update(count_rank='genus', count_label=genus, count_source='sequence_refined')
            elif genus and vcig >= o.sequence_fallback_vcig:
                row.update(count_rank='genus', count_label=genus, count_source='sequence_fallback')
            self.rows.append(row)
        changes = {}
        for row in self.rows:
            rank, label = row.get('count_rank'), row.get('count_label')
            if not rank:
                row['new_group'] = row['original_group']
                continue
            group = label if rank == 'species' else (f'GENUS_{label}' if rank == 'genus' else f'TAXON_{label}')
            row['new_group'] = group
            changes[row['contig']] = group
        return changes

    def record_strain_outcome(self, clusters):
        """clusters: StrainRefine's final strain clusters (lists of reference
        names). For contigs placed in a species: did the contig join a known
        strain's cluster, or form its own (a new strain)?"""
        where = {}
        for cl in clusters:
            for ref in cl:
                where[ref] = cl
        for row in self.rows:
            if row.get('count_rank') != 'species':
                continue
            cl = where.get(row['contig'])
            if cl is None:
                row['strain_outcome'] = 'no reads left on contig after strain clustering'
                continue
            known = [r for r in cl if r != row['contig'] and '|' in r]
            row['strain_outcome'] = 'joined known strain cluster (' + ','.join(known[:3]) + ')' if known else 'new strain (own cluster)'
            row['strain_cluster_representative'] = cl[0]

    COLUMNS = [
        'contig',
        'original_group',
        'new_group',
        'candidate_family',
        'best_hit',
        'best_hit_vcig',
        'contig_reads',
        'contig_final_reads',
        'resolved_rank',
        'resolved_label',
        'pass_used',
        'cluster_members',
        'contig_shared_fraction',
        'count_rank',
        'count_label',
        'count_source',
        'strain_outcome',
        'strain_cluster_representative',
        'decision',
    ]

    def write(self, out_tsv, final_reads, report_dir=None, names_path=None):
        for row in self.rows:
            row['contig_final_reads'] = final_reads.get(row['contig'], 0)
        with open(out_tsv, 'w') as fh:
            w = csv.DictWriter(fh, fieldnames=self.COLUMNS, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(self.rows)
        if not (report_dir and names_path):
            return
        o = self.opts
        report_dir = Path(report_dir)
        report_dir.mkdir(parents=True, exist_ok=True)
        names = load_names(names_path)
        neighbour_rows = []
        for row in self.rows:
            cid, fam = row['contig'], row['candidate_family']
            pools = self.species_pools(fam)
            best = row.get('best_hit') or ''
            neighbour_rows += write_contig_report(
                cid,
                fam,
                self.read_sets.get(cid, set()),
                pools,
                (row.get('count_rank') or '', row.get('count_label') or '', row.get('count_source') or ''),
                names,
                self.rank_of,
                o.cut_distance_pass2,
                o.cut_distance_pass3,
                o.min_population_size,
                report_dir / f'{cid}.png',
                best_hit=best.split('|')[1] if '|' in best else None,
                best_vcig=float(row.get('best_hit_vcig') or 0),
                pass_used=row.get('pass_used', ''),
                min_shared_reads=o.min_shared_reads,
                min_contig_fraction=o.min_contig_fraction,
                linkage_method=o.linkage_method,
                ref_pools={r: self.read_sets.get(r, set()) for r in self.family_refs.get(fam, [])},
                species_of=self.species_of,
                cut_refs=o.cut_distance_pass1,
                rules=dict(species_min_contig_fraction=o.species_min_contig_fraction, genus_min_vcig=o.genus_min_vcig,
                           sequence_fallback_vcig=o.sequence_fallback_vcig),
                shared_fraction=row.get('contig_shared_fraction'),
                resolved=(row.get('resolved_rank', ''), row.get('resolved_label', '')),
                final_reads=row['contig_final_reads'],
                single_genus_family=o.single_genus_family,
            )
        with open(report_dir / 'novel_contig_neighbours.tsv', 'w') as fh:
            cols = [
                'contig',
                'level',
                'neighbour_taxid',
                'neighbour_name',
                'neighbour_genus',
                'distance',
                'shared_reads',
                'contig_reads',
                'neighbour_reads',
                'fraction_of_contig_reads',
                'shown_in_figure',
            ]
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
            w.writeheader()
            w.writerows(neighbour_rows)
