#!/usr/bin/env python3
"""Prepare a reference database for MADRe.

MADRe needs, besides the taxonomy (nodes.dmp, names.dmp; the species of every reference comes from it at run time):
  1. a reference FASTA whose headers carry the taxid:  >|taxid|accession  (or >kraken:taxid|taxid|accession)
  2. a reference genome map (TSV): every sequence -> its genome, and that genome's total length

Two ways to build them:

  genomes   a folder with one FASTA per genome (e.g. NCBI downloads, own isolates, MAGs).
            Writes both files. Taxids from NCBI's assembly_summary.txt or from a genome -> taxid table.

  fasta     an existing combined FASTA that already has taxids in its headers (e.g. a Kraken library).
            Writes file 2; the FASTA is used as it is. Genomes from a sequence -> genome table,
            or grouped automatically (WGS project prefix for draft genomes, strain name for complete ones).

See DATABASE_PREPARATION.md for examples.
"""
import argparse
import csv
import gzip
import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

from NovelPlacement import ancestor, read_taxonomy  # the same taxonomy reading as MADRe itself

FASTA_SUFFIXES = ('.fna', '.fa', '.fasta', '.fas', '.fna.gz', '.fa.gz', '.fasta.gz', '.fas.gz')
ASSEMBLY_ACCESSION = re.compile(r'GC[AF]_\d+\.\d+')
WGS_ACCESSION = re.compile(r'^((?:NZ_)?[A-Z]{4,6}\d{2})\d{6,}\.\d+$')
REPLICON = re.compile(
    r'(,? (chromosome|plasmid|complete genome|complete sequence|genome assembly|whole genome shotgun|megaplasmid|contig|scaffold)\b.*$)', re.IGNORECASE
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.StreamHandler()])


def open_text(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def write_atomic(path, write):
    """Write to path + '.tmp' and rename when complete, so a failed run leaves no partial file."""
    tmp = f'{path}.tmp'
    write(tmp)
    os.replace(tmp, path)


def header_taxid(reference):
    """'|taxid|acc' or 'kraken:taxid|taxid|acc' -> taxid (the field StrainRefine reads)."""
    parts = reference.split('|')
    return parts[1] if len(parts) >= 3 and parts[1] else None


def fasta_records(path):
    """Yield (header line without '>', sequence lines) per record, streaming."""
    header, lines = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header, lines
                header, lines = line[1:].rstrip('\n'), []
            elif line.strip():
                lines.append(line if line.endswith('\n') else line + '\n')
    if header is not None:
        yield header, lines



def check_taxonomy(taxids, nodes):
    """Optional (--nodes): warn about database taxids missing from the taxonomy (MADRe maps them to themselves, i.e. no species)."""
    if not nodes:
        return
    parents, ranks = read_taxonomy(nodes)
    merged_path = Path(nodes).resolve().parent / 'merged.dmp'
    merged = {old: new for old, new, *_ in (line.replace('\t|', '').split() for line in open(merged_path))} if merged_path.exists() else {}
    n_merged = sum(1 for t in taxids if t not in parents and merged.get(t) in parents)
    if n_merged:
        logging.info(f'{n_merged} taxids were merged into newer taxids (merged.dmp); their species is taken from the newer taxid')
    taxids = {t if t in parents else merged.get(t, t) for t in taxids}
    missing = sorted(t for t in taxids if t not in parents)
    if missing:
        logging.warning(f'{len(missing)} taxids are not in {nodes}; MADRe will treat each as its own species, e.g. {", ".join(missing[:5])}')
    species = {ancestor(t, parents, ranks, 'species') or t for t in taxids}
    logging.info(f'Taxonomy check: {len(taxids)} taxids, {len(species)} species, {len(missing)} missing from {nodes}')


def write_genome_map(path, rows):
    """rows: [(reference, genome_id, sequence_length)]; genome_length = sum over the genome's sequences."""
    total = defaultdict(int)
    for _, genome, length in rows:
        total[genome] += length

    def write(tmp):
        with open(tmp, 'w', newline='') as handle:
            w = csv.writer(handle, delimiter='\t', lineterminator='\n')
            w.writerow(['reference', 'genome_id', 'genome_length'])
            for reference, genome, _ in rows:
                w.writerow([reference, genome, total[genome]])

    write_atomic(path, write)
    logging.info(f'Reference genome map: {len(rows)} sequences in {len(total)} genomes -> {path}')


def read_two_column(path, what):
    table = {}
    with open(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 2 or not f[0] or not f[1]:
                sys.exit(f'Error: {path} line {n}: expected two tab-separated columns ({what})')
            if f[0] in table and table[f[0]] != f[1]:
                sys.exit(f'Error: {path}: {f[0]} is given two different values ({table[f[0]]}, {f[1]})')
            table[f[0]] = f[1]
    return table


def read_assembly_summary(path):
    """NCBI assembly_summary.txt: assembly accession -> taxid (column 6; species taxid, column 7, if the taxid is missing)."""
    taxids = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 7:
                sys.exit(f'Error: {path}: assembly_summary.txt needs at least 7 columns')
            taxid = f[5] if f[5].isdigit() and int(f[5]) > 0 else f[6]
            if taxid.isdigit() and int(taxid) > 0:
                taxids[f[0]] = taxid
    return taxids


def genome_id_of(path):
    match = ASSEMBLY_ACCESSION.search(path.name)
    if match:
        return match.group()
    name = path.name
    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


# =============================================================================
# genomes: a folder with one FASTA per genome
# =============================================================================


def prepare_from_genomes(args):
    files = sorted(p for p in Path(args.genomes).rglob('*') if p.is_file() and p.name.endswith(FASTA_SUFFIXES))
    if not files:
        sys.exit(f'Error: no FASTA files ({", ".join(FASTA_SUFFIXES)}) found in {args.genomes}')
    taxid_of = read_assembly_summary(args.assembly_summary) if args.assembly_summary else read_two_column(args.genome_taxids, 'genome_id, taxid')
    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    fasta_out, map_out = Path(f'{out}.fna.gz'), Path(f'{out}.reference_genome_map.tsv')

    genome_ids, seen_sequences, rows, taxids, skipped = set(), set(), [], set(), []
    tmp = f'{fasta_out}.tmp'
    with gzip.open(tmp, 'wt', compresslevel=4) as handle:
        for i, path in enumerate(files, 1):
            genome = genome_id_of(path)
            if genome in genome_ids:
                sys.exit(f'Error: two FASTA files give the same genome id {genome} ({path.name}); use one file per genome')
            genome_ids.add(genome)
            taxid = taxid_of.get(genome)
            if taxid is None:
                skipped.append(genome)
                continue
            taxids.add(taxid)
            for header, lines in fasta_records(path):
                accession = header.split()[0] if header.split() else ''
                if not accession:
                    sys.exit(f'Error: empty FASTA header in {path}')
                accession = accession.replace('|', '_')
                if accession in seen_sequences:
                    sys.exit(f'Error: sequence {accession} occurs in more than one genome ({path.name}); sequence names must be unique')
                seen_sequences.add(accession)
                reference = f'|{taxid}|{accession}'
                handle.write(f'>{reference}\n')
                handle.writelines(lines)
                rows.append((reference, genome, sum(len(x.rstrip()) for x in lines)))
            if i == 1 or i % 1000 == 0 or i == len(files):
                logging.info(f'Processed {i}/{len(files)} genome files')
    if skipped:
        logging.warning(f'{len(skipped)} genomes have no taxid and were left out, e.g. {", ".join(skipped[:5])}')
    if not rows:
        os.remove(tmp)
        sys.exit('Error: no genome had a taxid; nothing written')
    os.replace(tmp, fasta_out)
    logging.info(f'Reference FASTA: {len(rows)} sequences from {len(genome_ids) - len(skipped)} genomes -> {fasta_out}')
    write_genome_map(map_out, rows)
    check_taxonomy(taxids, args.nodes)


# =============================================================================
# fasta: an existing combined FASTA with taxids in the headers
# =============================================================================


def auto_genome(reference, description):
    """Genome of a sequence when no sequence -> genome table is given:
    WGS draft accession -> taxid + WGS project prefix; otherwise taxid + strain name (header description without the
    replicon part, e.g. '... strain X chromosome' and '... strain X plasmid pY' -> '... strain X'); otherwise its own genome."""
    taxid, accession = header_taxid(reference), reference.split('|')[-1]
    match = WGS_ACCESSION.match(accession)
    if match:
        return f'{taxid}:{match.group(1)}', 'wgs project'
    strain = REPLICON.sub('', description).strip().rstrip(',')
    if strain:
        return f'{taxid}:{strain}', 'strain name'
    return f'{taxid}:{accession}', 'single sequence'


def prepare_from_fasta(args):
    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    map_out = Path(f'{out}.reference_genome_map.tsv')
    genome_of = read_two_column(args.sequence_genomes, 'sequence, genome_id') if args.sequence_genomes else None

    rows, taxids, seen, methods, no_taxid, no_genome = [], set(), set(), defaultdict(int), [], []
    for n, (header, lines) in enumerate(fasta_records(args.fasta), 1):
        reference, _, description = header.partition(' ')
        if reference in seen:
            sys.exit(f'Error: duplicate sequence name {reference} in {args.fasta}')
        seen.add(reference)
        taxid = header_taxid(reference)
        if taxid is None:
            no_taxid.append(reference)
            continue
        taxids.add(taxid)
        if genome_of is not None:
            genome = genome_of.get(reference) or genome_of.get(reference.split('|')[-1])
            if genome is None:
                no_genome.append(reference)
                continue
            methods['table'] += 1
        else:
            genome, method = auto_genome(reference, description)
            methods[method] += 1
        rows.append((reference, genome, sum(len(x.rstrip()) for x in lines)))
        if n % 100000 == 0:
            logging.info(f'Read {n} sequences')
    if no_taxid:
        sys.exit(
            f'Error: {len(no_taxid)} sequences have no taxid in the header (expected >|taxid|accession or >kraken:taxid|taxid|accession), '
            f'e.g. {", ".join(no_taxid[:3])}'
        )
    if no_genome:
        sys.exit(f'Error: {len(no_genome)} sequences are missing from {args.sequence_genomes}, e.g. {", ".join(no_genome[:3])}')
    logging.info('Genome assignment: ' + ', '.join(f'{k}: {v} sequences' for k, v in sorted(methods.items())))
    if methods.get('single sequence'):
        logging.warning(f'{methods["single sequence"]} sequences had neither a WGS accession nor a strain name and count as their own genome')
    write_genome_map(map_out, rows)
    check_taxonomy(taxids, args.nodes)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='mode', required=True)

    g = sub.add_parser('genomes', help='folder with one FASTA per genome -> reference FASTA, genome map')
    g.add_argument(
        '--genomes',
        required=True,
        help='Folder with one FASTA per genome (.fna/.fa/.fasta, optionally .gz), searched recursively. '
        'Genome id = NCBI assembly accession in the file name (GCF_/GCA_), otherwise the file name.',
    )
    source = g.add_mutually_exclusive_group(required=True)
    source.add_argument('--assembly-summary', help="NCBI assembly_summary.txt (taxid from column 6, or column 7 if missing)")
    source.add_argument('--genome-taxids', help='Tab-separated: genome_id <tab> taxid (genome_id as derived from the file name)')
    g.add_argument('--nodes', help='Optional: taxonomy nodes.dmp, to check that every database taxid is in it')
    g.add_argument('--out-prefix', required=True, help='Writes <prefix>.fna.gz and <prefix>.reference_genome_map.tsv')

    f = sub.add_parser('fasta', help='combined FASTA with taxids in headers -> genome map')
    f.add_argument('--fasta', required=True, help='Reference FASTA (optionally .gz) with headers >|taxid|accession or >kraken:taxid|taxid|accession')
    f.add_argument(
        '--sequence-genomes',
        help='Tab-separated: sequence <tab> genome_id (sequence = full header name or accession). '
        'Without it genomes are grouped automatically (WGS project prefix / strain name).',
    )
    f.add_argument('--nodes', help='Optional: taxonomy nodes.dmp, to check that every database taxid is in it')
    f.add_argument('--out-prefix', required=True, help='Writes <prefix>.reference_genome_map.tsv')

    args = parser.parse_args()
    prepare_from_genomes(args) if args.mode == 'genomes' else prepare_from_fasta(args)


if __name__ == '__main__':
    main()
