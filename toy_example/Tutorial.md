# How to run MADRe?

We prepared a short tutorial demonstrating how to run MADRe using a small example dataset. The reads (`sim_small.fastq.gz`) can be obtained from the [Zenodo page](https://zenodo.org/records/14845576). The example database, which contains only genomes from species included in the `sim_small` dataset, is available in this repository at `toy_example/db.fa.gz`.  
Download both of these files. This tutorial assumes that MADRe has been installed via Conda (see *Installation – OPTION 1* in main README file).

## Run the whole MADRe pipeline

The first step is to prepare the `config.ini` file:
```
[PATHS]
metaflye = path_to_flye
metaMDBG = path_to_metaMDBG
minimap = path_to_minimap2
hairsplitter = path_to_hairsplitter.py
seqkit = path_to_seqkit

[DATABASE]
predefined_db = MADRe/toy_example/db.fa.gz
strain_species_json = MADRe/database/taxids_species.json
```

For more information about MADRe’s pipeline parameters:
```
$ madre --help
usage: madre [-h] [--version] --out-folder OUT_FOLDER --reads READS [--reads_flag {pacbio,hifi,ont}] [--threads THREADS] [-F] [--config CONFIG] [--strictness {less-strict,strict,very-strict}] [--collapsed_strains_overhead COLLAPSED_STRAINS_OVERHEAD] [--min_contig_len MIN_CONTIG_LEN] [--use-myloasm]

MADRe

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --out-folder OUT_FOLDER
                        Path to the output folder.
  --reads READS         Path to the reads file (fastq/fq can be gzipped).
  --reads_flag {pacbio,hifi,ont}
                        Reads technology. (default=ont)
  --threads THREADS     Number of threads (default=32).
  -F, --force           Force rerun all steps.
  --config CONFIG       Path to the configuration file. (default=./config.ini)
  --strictness {less-strict,strict,very-strict}
                        Database reduction strictness level. (default=very-strict)
  --collapsed_strains_overhead COLLAPSED_STRAINS_OVERHEAD
                        Overhead for collapsed strains during database reduction. (default=2)
  --min_contig_len MIN_CONTIG_LEN
                        Minimum contig length for database reduction. (default=1000)
  --use-myloasm         Use Myloasm assembler instead of metaFlye/metaMDBG.

```

Once the `config.ini` file is ready, you can run the full pipeline with:

```
$ madre --out-folder MADRe_toy_example --reads sim_small.fastq.gz --reads_flag ont --threads 32 --config ./config.ini --strictness very-strict --collapsed_strains_overhead 2 --min_contig_len 1000
```

Detailed information about additional MADRe parameters can be found in the [Additional parameters information](#additional-parameters-information) section.

### Outputs
``` metaflye ``` - metaFlye output directory. For more information, check the [Flye GitHub page](https://github.com/mikolmogorov/Flye).

``` metaflye/hairsplitter ``` - HairSplitter output directory. For more information, check the [HairSplitter GitHub page](https://github.com/RolandFaure/Hairsplitter).
 
``` metaflye/collapsed_strains.txt ``` - Information about the number of collapsed strains per contig.

``` assembly.to_big_db.paf ``` - minimap2 output file from the first *Database Reduction* step. Contains mappings of assembled contigs to the initial database.

``` genomes_in_reduced.txt ``` - List of genome IDs selected for the reduced database.

``` reduced_db.fa ``` - Reduced database containing genomes listed in `genomes_in_reduced.txt`.

``` reads.to_reduced.paf ``` - minimap2 output file from the second *Read Classification* step. Contains mappings of reads to `reduced_db.fa`.

``` read_classification.out ``` - Output file from the *Read Classification* step. Each row represents the classification result for one read: `read_id : genome_id`.

``` rc_abundances.out ``` - Output file from the final *Calculate Abundances* step containing read count abundances. Each row represents the read count for a genome ID: `genome_id : read_count`.

``` abundances.out ``` - Output file from the final *Calculate Abundances* step containing abundances calculated as *sum_of_read_lengths/reference_length*. Each row represents abundance information for one genome ID: `genome_id : abundance`.

``` *.log ``` - Log files for various tools and steps.

### NOTE

All output files from this run could be found in ```MADRe/toy_example/madre_out```.

The default MADRe pipeline **does not** perform post-classification clustering. To perform clustering, follow the instructions under [Clustering](#clustering).  
**Please note that this feature is still under development!**

## Run *Database Reduction* step

*Database Reduction* step can be run independently, and the resulting reduced database can be used for other metagenomic analysis tasks.

As before, the first step is to prepare the `config.ini` file. Instructions for this can be found under [Run the whole MADRe pipeline](#run-the-whole-madre-pipeline).

For more information about MADRe’s *Database Reduction* parameters:

```
$ database-reduction --help
usage: database-reduction [-h] --database DATABASE --strain_species_info STRAIN_SPECIES_INFO --paf_path PAF_PATH --num_collapsed_strains NUM_COLLAPSED_STRAINS --reduced_list_txt REDUCED_LIST_TXT [--reduced_db REDUCED_DB] [--mapping_class MAPPING_CLASS] [--mapping_reduced_db MAPPING_REDUCED_DB] [--threads THREADS] [--strictness {less-strict,strict,very-strict}] [--min_contig_len MIN_CONTIG_LEN][--collapsed_strains_overhead COLLAPSED_STRAINS_OVERHEAD]

MADRe.

optional arguments:
  -h, --help            show this help message and exit
  --database DATABASE   Path to the strating database file (fasta/fna).
  --strain_species_info STRAIN_SPECIES_INFO
                        An additional parameter required if a custom database path is provided. JSON file with info about species taxid
                        for every strain taxid in the database. If you want to use default one provide path to
                        MADRe/database/taxids_species.json.
  --paf_path PAF_PATH   Path to the PAF file of assembly mapped to database.
  --num_collapsed_strains NUM_COLLAPSED_STRAINS
                        File containing info about number of collapsed strains for every contig (hairsplitter output).
  --reduced_list_txt REDUCED_LIST_TXT
                        Path to the file with list of genomes for reduced database.
  --reduced_db REDUCED_DB
                        Path to the reduced database file (fasta).
  --mapping_class MAPPING_CLASS
                        Path to the output mapping contig classification.
  --mapping_reduced_db MAPPING_REDUCED_DB
                        Path to the output mapping reduced database.
  --threads THREADS     Number of threads (default=32).
  --strictness {less-strict,strict,very-strict}
                        Strictness of database reduction - choices: less-strict, strict, very-strict - default: very-strict
  --min_contig_len MIN_CONTIG_LEN
                        Filter out contigs shorter than min_contig_len (default=1000).
  --collapsed_strains_overhead COLLAPSED_STRAINS_OVERHEAD
                        Maximum overhead for number of collapsed strains per contig estimated by HairSplitter (default=2).
```

To run this step, mappings of contigs to the database and the number of strains collapsed within each contig are required.  
These can be obtained by assembling the reads using, for example, `metaFlye`, and then running `minimap2` and `HairSplitter`:

```
$ flye --nano-raw sim_small.fastq.gz --out-dir metaflye_out --threads 32 --meta

$ minimap2 -x asm5 db.fa.gz metaflye_out/assembly.fasta -t 32 > assembly.to_big_db.paf

$ hairsplitter.py -i metaflye_out/assembly.fasta -f sim_small.fastq.gz -t 32 -o hairsplitter_out_dir -F
$grep '^>' hairsplitter_out_dir/tmp/cleaned_assembly.fasta | sed 's/^>//; s/@/:/; s/edge/contig/' > collapsed_strains.txt
```

Once `assembly.to_big_db.paf` and `collapsed_strains.txt` are generated, the *Database Reduction* step can be executed:

```
$ database-reduction --database db.fa.gz --strain_species_info MADRe/database/taxids_species.json --paf_path assembly.to_big_db.paf --num_collapsed_strains collapsed_strains.txt --reduced_list_txt reduced_list.txt --reduced_db reduced_db.fa --threads 32 --strictness very-strict --min_contig_len 1000 --collapsed_strains_overhead 2
```

Output files are described in the section [Run the whole MADRe pipeline – Outputs](#outputs).  
Additional parameters not used in the command above are explained in the section [Additional parameters information](#additional-parameters-information).

## Run *Read Classification* step

Just like the *Database Reduction* step, the *Read Classification* step can be run independently. For example, if the user already has prior knowledge about their sample and knows which organisms should be included in the database, they can build a smaller custom database and run only this step.  
In this tutorial, we will perform this step using the already reduced database available at `MADRe/toy_example/madre_out/reduced_db.fa`.

For more information about MADRe’s *Read Classification* parameters:

```
read-classification --help
read-classification [-h] --paf_path PAF_PATH --strain_species_info STRAIN_SPECIES_INFO [--read_class_output READ_CLASS_OUTPUT] [--clustering_out CLUSTERING_OUT]

MADRe.

optional arguments:
  -h, --help            show this help message and exit
  --paf_path PAF_PATH   Path to the PAF file of assembly mapped to database.
  --strain_species_info STRAIN_SPECIES_INFO
                        An additional parameter required if a custom database path is provided. JSON file with info about species taxid for every strain taxid in the database. If you want to use default one
                        provide path to MADRe/database/taxids_species.json.
  --read_class_output READ_CLASS_OUTPUT
                        Path to the output file with classification labels for reads. (default=read_classification.out)
  --clustering_out CLUSTERING_OUT
                        Path to clustering output directory. If provided clustering using mapping info from --paf_path will be performed, otherwise not.
```
To run this step, mappings of reads to the reduced database are required. These can be obtained using `minimap2`:

```
$ minimap2 -cx map-ont MADRe/toy_example/madre_out/reduced_db.fa sim_small.fastq.gz -t 32 > reads.to_reduced.paf
```
Once the `reads.to_reduced.paf` file is generated, we can run the *Read Classification* step:
```
read-classification --paf_path reads.to_reduced.paf --strain_species_info MADRe/database/taxids_species.json --read_class_output read_classification.out
```
Output files are described in the section [Run the whole MADRe pipeline – Outputs](#outputs).  
Additional parameters not used in the command above are explained in the section [Additional parameters information](#additional-parameters-information).

## Run *Calculate Abundances* step

The *Read Classification* step produces only an output file containing classification labels for each read.  
With the *Calculate Abundances* step, we can generate output files containing both read counts and relative abundances.

To run this step, the `read_classification.out` file from the *Read Classification* step is required (see [Run *Read Classification* step](#run-read-classification-step)).

For more information about MADRe’s *Calculate Abundances* parameters:
```
calculate-abundances --help
calculate-abundances [-h] [--db DB] --reads READS --read_class READ_CLASS [--rc_abundances_out RC_ABUNDANCES_OUT] [--abundances_out ABUNDANCES_OUT] [--clusters CLUSTERS]

MADRe.

optional arguments:
  -h, --help            show this help message and exit
  --db DB               Path to the database file (fasta). If DatabaseReduction used - path to the reduced database. WARNING: Required for estimated abundances calculation.
  --reads READS         Path to the reads file (fastq/fasta, can be gziped).
  --read_class READ_CLASS
                        Path to the input file with classification labels for reads from read classification step. (default=read_classification.out)
  --rc_abundances_out RC_ABUNDANCES_OUT
                        Path to the output file with read count. (default=rc_abundances.out)
  --abundances_out ABUNDANCES_OUT
                        Path to the output file with estimated abundances. If path is not given this file is not going to be generated. WARNING: In case of large sample and large database that can be computationally exhaustive job.
  --clusters CLUSTERS   Path to dir that contains clusters.txt and representatives.txt files. If provided, the abundances in output files will be reported including cluster information.
```
The *Calculate Abundances* step can be run as:
```
$ calculate-abundances --reads sim_small.fastq.gz --read_class read_classification.out --rc_abundances_out rc_abundances.out
```

Relative abundances are calculated as: `sum_of_read_lengths / reference_length`.  
This calculation can be computationally demanding, so it is **not** included in the default *Calculate Abundances* run (**but it is part of the full MADRe pipeline**).  
If you also want to generate the output file with relative abundances, run:

```
$ calculate-abundances --reads sim_small.fastq.gz --read_class read_classification.out --rc_abundances_out rc_abundances.out --db db.fa.gz --abundances_out abundances.out
```
Output files are described in the section [Run the whole MADRe pipeline – Outputs](#outputs).  
Additional parameters not used in the command above are explained in the section [Additional parameters information](#additional-parameters-information).

## Clustering

The current implementation of MADRe’s independent steps can also be used to obtain clustering of very similar strains based on their mapping profiles, as utilized in the *Read Classification* step.

Please note that this feature is still under development and is not part of default MADRe's pipeline.

Cluster data can be obtained by running the *Read Classification* step as follows:

```
$ read-classification --paf_path reads.to_reduced.paf --strain_species_info MADRe/database/taxids_species.json --read_class_output read_classification.out --clustering_out clusters 
```
`clusters` is the output directory containing two files: `clusters.txt` and `representatives.txt`.  
Example outputs can be found in `MADRe/toy_example/madre_out/clusters`.  

Each line in `clusters.txt` represents a single cluster — all genome IDs belonging to that cluster are separated by spaces.  
The file `representatives.txt` contains the representative genomes for each cluster.  
The lines in these two files correspond to the same clusters, e.g., the representative genome for the cluster listed in line 7 of `clusters.txt` will be found in line 7 of `representatives.txt`.

After this step, the `read_classification.out` file will remain unchanged and will contain the same labels as when clustering was not performed.

However, by combining the files from the `clusters` directory with the *Calculate Abundances* step, we can compute cluster-level abundances:

```
$ calculate-abundances --reads sim_small.fastq.gz --read_class read_classification.out --rc_abundances_out rc_abundances.with_clusters.out --db db.fa.gz --abundances_out abundances.with_clusters.out --clusters clusters
```

After running this command, the output files `rc_abundances.with_clusters.out` and `abundances.with_clusters.out` will contain entries in the format  `genome_ID_with_cluster : abundance`, representing the abundance of each cluster, with the cluster identified by the representative genome listed with that genome ID.  
Example output files can be found in `MADRe/toy_example/madre_out`.

## Additional parameters information 

Default parameter values produce successful results on our benchmarking datasets. However, real metagenomic samples are often highly complex, and proper metagenome analysis may sometimes require parameter tuning.

This section provides an explanation of several additional MADRe parameters and their effects:

- `--strictness` 

This parameter can be set when running `madre` or `database-reduction`. It defines how strict the database reduction process will be.  
There are three options: `very-strict`, `strict`, and `less-strict`. The default option is `very-strict`. Using less strict options may increase the number of false-positive identifications in the reduced database but can also help include low-abundance strains that might otherwise be excluded under stricter settings.

- `--collapsed_strains_overhead`

This parameter can be set when running `madre` or `database-reduction`. It defines how many more strains (beyond the number estimated by HairSplitter for a given contig) can be considered when evaluating mapping scores for database reduction.  
The default value is `2`. Increasing this value has a similar effect to reducing the strictness, potentially including more genomes but also increasing the chance of false positives.

- `--min_contig_len` 

This parameter can be set when running `madre` or `database-reduction`. All contigs shorter than the specified value are filtered out.  
The default value is `1000`. Very small contigs with unusually high coverage often represent misassemblies and can introduce noise into the database reduction process. 

- `--use-myloasm`

This parameter can be set when running `madre` or `database-reduction`. The default value is *False*.  
If set to *True*, MADRe will use the *Myloasm* assembler instead of *metaFlye* or *metaMDBG*.  
Our benchmarking results showed comparable performance, with *metaFlye* performing better for high-abundance strains and *Myloasm* for low-abundance strains.  
Since *“there is no free lunch,”* different assembly approaches can help address different challenges. For this reason, we plan to include additional assemblers in future versions of the MADRe pipeline.

- `--strain_species_info`

This parameter can be set when running `madre`, `database-reduction`, or `read-classification`.  
It specifies the path to a file such as `MADRe/database/taxids_species.json`, which contains information mapping each strain, subspecies, and species taxID to its corresponding species-level taxID.  
If database contains taxID that is not part of `MADRe/database/taxids_species.json` file, MADRe pipeline will continue but it will skip mappings containing that organism which can impact classification results.  
The default JSON file contains taxonomy information from December 2024.  
To update to a newer taxonomy version, a new `taxids_species.json` file can be generated using the script `MADRe/database/metagenomika_eval/build_taxids_index.py`.

- `--mapping-class` 

This parameter can be set when running `database-reduction`.  
The corresponding output file will contain information on which reference genome each contig has the highest mapping score against.

- `--mapping-reduced-db` 

This parameter can be set when running `database-reduction`.  
The corresponding output file will contain the reduced database, including only reference genomes that the contigs mapped to with the highest scores.




