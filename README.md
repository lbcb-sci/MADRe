# MADRe

<p align="center">
  <img src="imgs/MADRe_logo-modified.png" alt="MADRe logo" width="150"/>
</p>


Strain-level metagenomic classification with Metagenome Assembly driven Database Reduction approach

+ [Why MADRe?](#why-madre)
+ [Installation](#installation)
+ [Build database](#build-database)
+ [How to run MADRe?](#how-to-run-madre)
+ [Run specific steps](#run-specific-steps)
+ [Abundance calculation](#abundance-calculation)
+ [Citing MADRe](#citing-madre)

## Why MADRe?

**MADRe (Metagenomic Assembly-Driven Database Reduction)** is designed for metagenomic analyses where there is **no prior knowledge about the sample composition** and the starting database is **large and diverse**, containing thousands of species and strains.

In such exploratory settings, traditional read-based classifiers either require extensive computational resources or struggle to resolve closely related genomes.  
MADRe overcomes these limitations by introducing an **assembly-guided database reduction strategy** that automatically identifies and retains only the genomes supported by the data, thereby enabling a more computationally efficient mapping-based classification process.  
This dramatically reduces both **runtime and disk usage** compared to traditional mapping-based classifiers, while improving **classification precision and accuracy** relative to k-mer-based metagenomic classification methods.

### When to use MADRe?
Use MADRe when working with:
- Complex metagenomic datasets where the taxonomic composition is unknown.
- Very large reference databases containing multiple strains per species.
- Long-read sequencing data (ONT, PacBio HiFi) where assembly is feasible.

### Why MADRe is different?
- **Efficient exploration of large databases** – Instead of mapping every read to every genome, MADRe narrows the search space through an assembly-driven reduction step, lowering computational load without significantly sacrificing accuracy. 
- **Resource-aware design** – For smaller datasets (~1.7 M ONT reads), MADRe requires up to ~2.5× less RAM and achieves ~5.2× shorter runtime, while for larger datasets (~5 M ONT reads) it runs up to ~3× faster and uses ~7.5× less disk space, all while maintaining higher interpretability and accuracy compared with other mapping-based, strain-aware classifiers.
- **Improved precision over k-mer based tools** – By leveraging alignment-based evidence from assembled contigs, MADRe avoids many of the false-positive assignments typical for k-mer classifiers.  
- **Modular and transparent** – Each step (*Database Reduction*, *Read Classification*, *Calculate Abundances*) can be executed independently, producing interpretable outputs suitable for downstream analyses.

MADRe is particularly useful as a **first-pass classification tool** for large, uncharacterized metagenomic datasets, providing a computationally efficient and biologically meaningful starting point for deeper strain-level analysis.

## Installation

### OPTION 1 : Conda

```
conda install bioconda::madre
```
set up the configuration (config.ini file):
```                                                               
[PATHS]
metaflye = flye
metaMDBG = metaMDBG
minimap = minimap2
hairsplitter = hairsplitter.py
seqkit = seqkit

[DATABASE]
predefined_db = /path/to/database.fna
strain_species_json = /path/to/taxids_species.json
```
NOTE: Prebuilt version of ```taxids_species.json``` can be found in GitHub database folder. More information about it find under the section [Build database](#build-database).

simple run:

```
madre --reads [path_to_the_reads] --out-folder [path_to_the_out_folder] --config config.ini
```
more information:
```
madre --help
```

### OPTION 2: Running from source

```
git clone https://github.com/lbcb-sci/MADRe
cd MADRe
```

For running from source you need to install following dependecies:
- python >= 3.10
- scikit-learn
- minimap2
- flye
- metamdbg
- hairsplitter
- seqkit
- kraken2
- myloasm (optional)

Dependencies can be installed through conda:
```
conda create -n MADRe_env python=3.10 scikit-learn minimap2 flye metamdbg hairsplitter seqkit kraken2 -c conda-forge -c bioconda 
conda activate MADRe_env
```

set up the configuration (config.ini file):
```                                                               
[PATHS]
metaflye = /path/to/flye
metaMDBG = /path/to/metaMDBG
minimap = /path/to/minimap2
hairsplitter = /path/to/hairsplitter.py
seqkit = /path/to/seqkit
myloasm = /path/to/myloasm #optional

[DATABASE]
predefined_db = /path/to/database.fna
strain_species_json = ./database/taxids_species.json
```

simple run:
```
python MADRe.py --reads [path_to_the_reads] --out-folder [path_to_the_out_folder] --config config.ini
```

more information:
```
python MADRe.py --help
```

The recommended database is Kraken2 bacteria database - instructions on how to build it you can find under the section [Build database](#build-database).

Information on how to run specific MADRe steps find under the section [Run specific steps](#run-specific-steps).
</details>

## Build database

### Recommended database (kraken2 built database)
The recommend database is the kraken2 built bacteria database following next steps:
```
kraken2-build --download-taxonomy --db $DBNAME
kraken2-build --download-library bacteria --db $DBNAME
kraken2-build --build --db $DBNAME
```

Detailed instructions that are including the one listed here can be found at [kraken2 github page](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown).

### GTDB database

For using GTDB database, first download the latest GTDB database version and its associated metadata from https://data.gtdb.aau.ecogenomic.org:
```
wget https://data.gtdb.aau.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
wget https://data.gtdb.aau.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
gunzip bac120_metadata.tsv.gz
```
Then run script `database/gtdb_to_madre.sh`:
```
./gtdb_to_madre.sh --tar gtdb_genomes_reps.tar.gz --meta bac120_metadata.tsv --out MADRe_reference_database
```

### Build your own database

If you want to use your database it is important to have taxonomy information for the references included in the database. 

References in the database should have headers in this way:

```
>|taxid|accession_number
```

```../database/taxids_species.json``` file contains information on species taxid for every strain taxid obtained from NCBI taxonomy (downloaded December 2024.). 

MADRe for species-level classification step uses taxids index. For building new taxids index from newer taxonomy or for different taxonomic levels you will need taxonomy files (can be downloaded here https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) and you can use ```database/build_json_taxids.py``` script.

## How to run MADRe?

This README contains basic information on how to run MADRe pipeline. However, **for a more detailed tuorial check [toy_example/Tutorial.md](https://github.com/lbcb-sci/MADRe/blob/main/toy_example/Tutorial.md) file**.

## Run specific steps

MADRe is the pipeline contained of two main steps: 1) database reduction and 2) read classification.

It is possible to run those steps independently. More infromation on running can be obtained with:

```
database-reduction --help
read-classification --help
```

installed from source:
```
python src/DatabaseReduction.py --help
python src/ReadClassification.py --help
```

### Database reduction information

To run database reduction step separately you need to provide names of the output paths, mapping PAF file containg contigs mappings to large database (database needs to follow rules from [Build database](#build-database) section) and text file containing how many strains are collapsed in which contig. If contig represents only one strain there should be 0 next to it, if it represents 2 strains, 1 is collapsed so there should be 1 next to it. The file should look like this:
```
...
contig_7:0 
contig_8:0 
contig_8:1 
contig_8:2 
contig_8:3
...
```

If as output you only specify ```--reduced_list_txt``` you won't get fasta file of reduced database, just list of references that should go to reduced database. To get fasta file of reduced database specify ```--reduced_db```.

Database reduction step uses taxid index. By default it uses ```database/taxid_species.json```. If specific large database is used, then right taxid index should be provided using ```--strain_species_info```.


### Read classification information

To run read classification step separately you need to provide PAF file containing read mappings to the reference. This step can be run on any database (database needs to follow rules from [Build database](#build-database) section), so it doesn't have to be previously reduced.

Read classification step uses taxid index. By default it uses ```database/taxid_species.json```. If specific large database is used, then right taxid index should be provided using ```--strain_species_info```.

Output file is text file containg lines as: ```read_id : reference```.

### Read Classification with clustering

As part of read classification step, clustering of very similar strains can also be performed. If you want to perform clustering provide path to the directory with output clustering files using ```--clustering_out```. Output clustering files are:
```
clusters.txt - Every line represents one cluster. References in cluster separated with spaces.
representatives.txt - Every line represents a cluster representative reference of the cluster from that line in clusters.txt file.
```

## Abundance calculation

For abundance calculation information run:
```
calculate-abundances --help
```
installed from source:
```
python src/CalculateAbundances.py --help
```
The input to this step is read classification output file that has lines as ```read_id : reference```. This file can be obtained with read classification step.

The default output is rc_abundances.out containing read count abundances. If you want to calculate abundance as sum_of_read_lengths/reference_length you need to provide database path used in read classification step using ```--db``` - be aware that this step if database is big takes a little bit longer than calculation of just read count abundances. 

If you want to calculate cluster abundances, you need to provide path to the directory containing ```clusters.txt``` and ```representatives.txt``` files. In that case output files will contain only represetative references with sumarized abundances for cluster that reference is represetative of.

## Citing MADRe
bioRxiv preprint - https://www.biorxiv.org/content/10.1101/2025.05.12.653324:
```
Lipovac, J., Sikic, M., Vicedomini, R., & Krizanovic, K. (2025). MADRe: Strain-Level Metagenomic Classification Through Assembly-Driven Database Reduction. bioRxiv, 2025-05.
```

