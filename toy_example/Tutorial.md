# MADRe toy example

This tutorial runs MADRe on a small simulated ONT dataset twice:

1. with a database that contains every species of the sample, and
2. with a database from which one species was removed, to show how MADRe reports an organism that is missing from the database.

In each run the assembly with metaFlye takes a few minutes; all other steps together take well under a minute. The tutorial assumes MADRe is installed (see the main [README](../README.md#installation)) and that you are in the `toy_example` folder.

## The data

**Reads:** `sim_small.fastq.gz` (6,358 simulated ONT reads) from [Zenodo](https://zenodo.org/records/14845576):

```bash
wget -O sim_small.fastq.gz "https://zenodo.org/records/14845576/files/sim_small.fastq.gz?download=1"
```

The reads come from four strains of two species, plus 130 random and junk reads:

| species | strain (taxid) | reads |
|---|---|---|
| *Streptococcus anginosus* | subsp. *anginosus* FDAARGOS_1569 (1671923) | 2,534 |
| | C1051 (862970) | 512 |
| *Adlercreutzia equolifaciens* | subsp. *celatus* JCM 14811 (394340) | 2,378 |
| | DSM 19450 (1384484) | 804 |

**Databases** (in this folder):

| file | content |
|---|---|
| `db.fa.gz` | 21 genomes: 19 *S. anginosus* strains and both *A. equolifaciens* strains |
| `db_without_adlercreutzia.fa.gz` | the same without *A. equolifaciens* (19 genomes) |

Each database comes with its genome map (`*.reference_genome_map.tsv`), made with `madre-prepare-database fasta --fasta db.fa.gz --out-prefix db`.

**Taxonomy:** `taxonomy/nodes.dmp` and `taxonomy/names.dmp` are a small extract of the NCBI taxonomy with only the lineages of these genomes, so the example needs no taxonomy download. For real data use the full NCBI taxonomy (see [DATABASE_PREPARATION.md](../DATABASE_PREPARATION.md)).

## config.ini

Write one config per database, with the paths of your tools (absolute paths):

```ini
[PATHS]
metaflye = flye
metaMDBG = metaMDBG
minimap = minimap2
hairsplitter = hairsplitter.py

[DATABASE]
predefined_db = /path/to/MADRe/toy_example/db.fa.gz
```

Save it as `config_full.ini`, and a copy with `predefined_db = /path/to/MADRe/toy_example/db_without_adlercreutzia.fa.gz` as `config_without.ini`.

## Run 1: all species in the database

```bash
madre --out-folder toy_full --reads sim_small.fastq.gz --reads_flag ont --threads 64 \
      --config config_full.ini --reference_genome_map db.reference_genome_map.tsv \
      --nodes taxonomy/nodes.dmp --names taxonomy/names.dmp
```

All 75 contigs are KNOWN (`toy_full/contig_classification.tsv`): every contig is well covered by database genomes.

**Species** (`toy_full/strainrefine_taxonomic_profile.tsv`, rank `species`):

| species | reads | genome copies | relative abundance |
|---|---|---|---|
| *Streptococcus anginosus* | 3,037 | 24.4 | 59.7% |
| *Adlercreutzia equolifaciens* | 3,182 | 16.5 | 40.3% |

The relative abundance is computed from genome copies, not from read counts: the *A. equolifaciens* genomes are longer, so the same number of reads means fewer copies.

**Strains** (`toy_full/strainrefine_reference_summary.tsv`):

| strain | true reads | assigned reads |
|---|---|---|
| *S. anginosus* subsp. *anginosus* FDAARGOS_1569 | 2,534 | 2,528 |
| *A. equolifaciens* subsp. *celatus* JCM 14811 | 2,378 | 2,378 |
| *A. equolifaciens* DSM 19450 | 804 | 804 |
| *S. anginosus* C1051 | 512 | 509 |

Of the 19 *S. anginosus* genomes in the database, the two that are in the sample receive the reads. `toy_full/strainrefine_clustering/clusters.txt` shows that FDAARGOS_1569 was clustered with two other *S. anginosus* genomes that the reads cannot tell apart.

## Run 2: *A. equolifaciens* missing from the database

```bash
madre --out-folder toy_without --reads sim_small.fastq.gz --reads_flag ont --threads 64 \
      --config config_without.ini --reference_genome_map db_without_adlercreutzia.reference_genome_map.tsv \
      --nodes taxonomy/nodes.dmp --names taxonomy/names.dmp --placement-report
```

`--placement-report` adds the placement figures explained below; they are not drawn by default.

**Contigs** (`toy_without/contig_classification.tsv`): 21 contigs are KNOWN (*S. anginosus*) and 53 are NOVEL_CANDIDATE. These are the *A. equolifaciens* contigs, which no database genome covers (`vcig_norm` close to 0).

**Novel groups** (`toy_without/novel_groups.tsv`): the 53 novel contigs are grouped into putative organisms by sequence similarity and assembly-graph links. Most of the genome forms one group, `NOVEL_COMPARE_00002`; short contigs that are not linked to it stay as small groups of their own.

**Reads:** 3,166 of the 3,182 *A. equolifaciens* reads are assigned to the novel contigs (`toy_without/read_classification.out`), instead of being forced onto *S. anginosus*. The *S. anginosus* result hardly changes: 3,032 reads, 50.4%, with the same two strains.

**Profile** (`toy_without/strainrefine_taxonomic_profile.tsv`):

| taxon or group | rank | reads | relative abundance |
|---|---|---|---|
| *Streptococcus anginosus* | species | 3,032 | 50.4% |
| `NOVEL_COMPARE_00002` | novel_group | 2,951 | 46.0% |
| 21 smaller novel groups | novel_group | 230 | 3.6% |

**Placement** (`toy_without/novel_placement.tsv`): MADRe tries to place each novel contig in the taxonomy by comparing it with its nearest database genomes. Here nothing related to *A. equolifaciens* is in the database: three contigs had weak hits to *Streptococcus* genomes, which suggested the family *Streptococcaceae*, but these hits did not hold up (`no_match_at_any_level`), so the contigs stay unplaced. The novel organism is therefore reported as a novel group rather than under a wrong taxon. With a relative of the missing organism in the database (e.g. another species of the same genus), the group would be placed at that genus.

**Placement figures** (`toy_without/novel_placement_report/contig_12.png`, `contig_45.png`, `contig_8.png`): one figure per contig with a candidate family. In `contig_45.png`:

- the header says "1224 reads → NOT COUNTED at any rank": the contig's reads are not counted under any taxon, so they stay in the novel group;
- the top panel (reference view) shows the contig next to the eight nearest *S. anginosus* genomes, clustered by read sharing (distance 1 − shared reads / reads of the smaller set). The contig shares at most 9 of its 1,318 reads with any of them and stays alone at the cut (red dashed line), so it gets no species;
- the two lower panels repeat the comparison against the species and the genus. With only *S. anginosus* and *Streptococcus* in the database there is nothing to compare at the species level (pass skipped), and the contig again stays alone at the genus level.

The README shows a figure of a contig that was placed, in genus *Helicobacter* ([Placement figures](../README.md#placement-figures)).

## Everything in the output folder

The tables above are MADRe's results. The other files are what each step passes on to the next; they are kept so that an interrupted run can be resumed (finished steps are skipped unless `-F` is given) and so that every decision can be traced. Listed in the order of the pipeline, for `toy_without/`:

**Results**

| file | content |
|---|---|
| `strainrefine_taxonomic_profile.tsv` | the taxonomic profile: one row per taxon at every rank, and per novel group; read count, genome copies (`genome_equivalent_signal`), relative abundance (%). With `--min-abundance`, the complete profile is in `strainrefine_taxonomic_profile.all_species.tsv`. |
| `strainrefine_species_or_group_summary.tsv` | one row per species or novel group: reads, bases, genome copies, relative abundance and the references or contigs behind it |
| `strainrefine_reference_summary.tsv` | one row per strain (genome) or contig with reads: read count, bases, length, depth, horizontal coverage; for a contig also the genome of its best database match (`best_database_hit`) and the contig's novelty score (`best_database_hit_vcig`, the `vcig_norm` of `contig_classification.tsv`); a label only |
| `read_classification.out` | `read_id : reference` for every assigned read: a database sequence, a contig of a novel organism, or the host. Reads not listed are unassigned. |

**Step 1, assembly** (`metaflye/`, or `metaMDBG/` / `myloasm_*/`)

| file | content |
|---|---|
| `metaflye/assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt` | the contigs, the assembly graph and per-contig coverage and circularity (the graph links are used to group novel contigs) |
| `metaflye/hairsplitter/`, `metaflye/collapsed_strains.txt` | step 3: HairSplitter's output, and the number of strains collapsed into each contig (`contig:number`) |
| other files in `metaflye/` | metaFlye's own intermediate folders and log |

**Steps 2–4, database reduction and contig classification**

| file | content |
|---|---|
| `assembly.to_big_db.paf` | the contigs mapped to the full database (minimap2) |
| `minimap_reduction.log` | minimap2's log of that mapping |
| `mapping_contig_class.txt` | `contig : reference`: the best-scoring database sequence of every contig with a database hit |
| `contig_classification.tsv` | every contig: KNOWN / UNCERTAIN / NOVEL_CANDIDATE, its novelty score (`vcig_norm`: the fraction of the contig covered by alignments to database genomes, weighted by their identity; 1 = fully known, 0 = nothing similar), best database hit and length |
| `reduced_db_refs.txt` | the database sequences selected for the reduced database (whole genomes) |
| `strain_species.json` | the species of every database taxid, read from `--nodes` |

**Step 5, novel organisms**

| file | content |
|---|---|
| `extra_sequences.fa` | the UNCERTAIN and NOVEL_CANDIDATE contigs |
| `extra_sequences.self.paf` | these contigs mapped against each other, for grouping |
| `novel_groups.tsv` | every novel contig with its groups: the strict group (`NOVEL_…`, near-identical contigs), the comparison group (`NOVEL_COMPARE_…`, contigs of one putative organism), how it was grouped, coverage, circularity and best database hit |
| `novel_id_map.json`, `novel_comparison_group_map.json` | the same two groupings as contig → group maps |
| `extended_strain_species.json` | `strain_species.json` extended with the contigs: UNCERTAIN contigs → the species of their best hit, novel contigs → their comparison group |

**Steps 6–7, extended database and read mapping**

| file | content |
|---|---|
| `extended_db.fa` | the reduced database + the UNCERTAIN and novel contigs (+ the host genome with `--host-reference`) |
| `extended_db.fa.inputs` | what `extended_db.fa` was built from (database, checksums of the selected references and contigs, host); a rerun rebuilds it only if these change |
| `reads.to_extended.paf` | the reads mapped to the extended database (minimap2) |
| `minimap_classification.log` | minimap2's log of that mapping |
| `reduced_db.fa` | only with `--keep-reduced-db`: the selected database sequences on their own |

**Steps 8–9, StrainRefine and abundances**

| file | content |
|---|---|
| `novel_placement.tsv` | the placement of every novel contig with a candidate family: the rank and taxon reached, which comparison decided it, or why it stayed unplaced (e.g. `no_match_at_any_level`) |
| `novel_placement_report/contig_*.png`, `novel_contig_neighbours.tsv` | only with `--placement-report`: one figure per such contig with its nearest database genomes, and the numbers behind the figures |
| `extended_strain_species.placed.json` | `extended_strain_species.json` after placement: placed contigs now point to their species, genus or family |
| `strainrefine_clustering/clusters.txt`, `representatives.txt` | groups of strains the reads cannot tell apart (one group per line) and the strain reported for each |

## Expected output

`expected_output/full_db/` and `expected_output/without_adlercreutzia/` hold the result tables of both runs (the read assignments compressed). Assembly with metaFlye is not fully deterministic with many threads, so your contigs, and with them the numbers, can differ slightly.
