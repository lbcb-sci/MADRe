# Preparing a database for MADRe

MADRe (`MADRe.py`) needs three inputs besides the reads:

| input | what it is | given as |
|---|---|---|
| reference FASTA | all reference genomes; every header carries the taxid | `predefined_db` in `config.ini` |
| reference genome map | every sequence → its genome, and that genome's total length | `--reference_genome_map` |
| taxonomy | `nodes.dmp` and `names.dmp` (NCBI format); the species of every reference is taken from it | `--nodes`, `--names` |

`PrepareDatabase.py` builds the first two. The taxonomy is downloaded from NCBI:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxdump && tar -xzf taxdump.tar.gz -C taxdump      # gives taxdump/nodes.dmp and taxdump/names.dmp
```

Use the taxonomy version the database's taxids come from.

## Which mode to use

- **You have a folder with one FASTA file per genome** (NCBI downloads, your own isolates, MAGs) → [mode `genomes`](#mode-genomes-one-fasta-per-genome). Builds both files.
- **You already have one combined FASTA with taxids in the headers** (e.g. a Kraken2 library) → [mode `fasta`](#mode-fasta-an-existing-combined-fasta). Builds the genome map; the FASTA is used as it is.

## Mode `genomes`: one FASTA per genome

```bash
python PrepareDatabase.py genomes \
    --genomes my_genomes/ \
    --assembly-summary assembly_summary.txt \
    --nodes taxdump/nodes.dmp \
    --out-prefix db/my_db
```

Writes `db/my_db.fna.gz` and `db/my_db.reference_genome_map.tsv`.

- **Genome files:** `.fna`, `.fa`, `.fasta` or `.fas`, optionally gzipped, searched recursively. Each file is one genome (chromosomes, plasmids, contigs of a draft genome).
- **Genome id:** the NCBI assembly accession in the file name if there is one (`GCF_000005845.2_ASM584v2_genomic.fna.gz` → `GCF_000005845.2`), otherwise the file name without its extension (`isolateA.fasta` → `isolateA`).
- **Taxids**, one of:
  - `--assembly-summary`: NCBI's `assembly_summary.txt` (e.g. `https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt`); the taxid is column 6, or the species taxid (column 7) when column 6 is empty;
  - `--genome-taxids`: a tab-separated table for genomes that are not in NCBI's summary:

    ```
    isolateA	562
    isolateB	28901
    ```

  Genomes without a taxid are left out, with a warning.
- **Sequence names** must be unique across all genomes. Headers are rewritten to `>|taxid|accession` (accession = the first word of the original header).

## Mode `fasta`: an existing combined FASTA

```bash
python PrepareDatabase.py fasta \
    --fasta library.fna \
    --nodes taxdump/nodes.dmp \
    --out-prefix db/library
```

Writes `db/library.reference_genome_map.tsv`. Use `library.fna` itself as `predefined_db`.

- **Headers** must carry the taxid as the second `|`-separated field: `>|562|NC_000913.3 ...` or `>kraken:taxid|562|NC_000913.3 ...`.
- **Genomes**, one of:
  - `--sequence-genomes`: a tab-separated table sequence → genome (sequence = the full header name or just the accession):

    ```
    kraken:taxid|562|chromosome_of_strainA	strainA
    kraken:taxid|562|plasmid_of_strainA	strainA
    kraken:taxid|562|chromosome_of_strainB	strainB
    ```

  - automatic, when no table is given:
    - draft genomes: sequences of the same WGS project (`NZ_ABCD01000123.1` → `NZ_ABCD01`) and taxid form one genome;
    - complete genomes: sequences with the same strain name in the header and the same taxid form one genome ("*E. coli* strain K-12 chromosome" and "*E. coli* strain K-12 plasmid F" → "*E. coli* strain K-12");
    - a sequence with neither is its own genome (counted in the log).

  Use the table whenever you have it; the automatic grouping relies on header descriptions.

## The output files

**Reference genome map** (`--reference_genome_map`): tab-separated, a header line, one row per sequence. `genome_length` is the whole genome's length (all its sequences), repeated on each of its rows. `reference` matches the FASTA header exactly (up to the first space).

```
reference	genome_id	genome_length
|1000373|NC_016757.1	GCF_000895895.1	17078
|1000373|NC_016758.1	GCF_000895895.1	17078
|562|NZ_CM137271.1	GCF_049197725.1	5411416
```

MADRe uses it to put whole genomes into the reduced database, to treat the sequences of one genome as one strain, and to divide reads by the whole genome's length.

## Large databases: index for fast extraction

For every sample, MADRe extracts the references it selected from the database. With a plain or gzip-compressed FASTA this reads the whole file each time (`seqkit grep`), which takes hours for a database of hundreds of GB. Compress it with `bgzip` and index it once; MADRe then reads the selected references directly (`samtools faidx`), in minutes:

```bash
zcat my_db.fna.gz | bgzip -@ 32 -c > my_db.bgz.fna.gz
samtools faidx my_db.bgz.fna.gz      # writes my_db.bgz.fna.gz.fai and .gzi next to it
```

Use the bgzip file as `predefined_db`. The extracted sequences are the same; only header descriptions after the first word are not copied. MADRe falls back to `seqkit grep` when there is no index, or when a sequence name ends in something that looks like a region (`name:100` or `name:100-200`).

The index also speeds up the mapping of the contigs to the database: minimap2 indexes a large database in parts of ~8 Gbp, one after another, and MADRe then maps several parts at a time (`--db-map-jobs`, default 8), with the same result as one minimap2 run. On a 653 Gbp RefSeq database this took 46 min instead of about 3 h.

## config.ini

```ini
[PATHS]
metaflye = flye
metaMDBG = metaMDBG
minimap = minimap2
hairsplitter = /path/to/Hairsplitter/hairsplitter.py
myloasm = myloasm

[DATABASE]
predefined_db = /path/to/db/my_db.fna.gz
```

Use absolute paths.

## Running MADRe

```bash
python MADRe.py \
    --out-folder results/sample1 \
    --reads sample1.fastq.gz \
    --reads_flag ont \
    --threads 32 \
    --config config.ini \
    --reference_genome_map /path/to/db/my_db.reference_genome_map.tsv \
    --nodes taxdump/nodes.dmp \
    --names taxdump/names.dmp
```

Optional: `--host-reference human.fna --host-name human` to separate host reads (e.g. human T2T-CHM13).

Optional: `--min-abundance X` leaves species below X% out of `strainrefine_taxonomic_profile.tsv`; the complete profile is kept as
`strainrefine_taxonomic_profile.all_species.tsv`. Off by default.

## Checks done while preparing

`PrepareDatabase.py` stops with a message if a header has no taxid, a sequence name occurs twice, or a sequence is missing from `--sequence-genomes`. It warns about genomes without a taxid (left out) and taxids missing from `nodes.dmp` (mapped to themselves). Output files are written under a temporary name and renamed when complete.
