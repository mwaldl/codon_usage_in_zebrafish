# analysing codon usage in zebrafish

pipeline used to analyse codon usage during zebrafish embryo development as presented in "tRNA expression and modification landscapes, and their dynamics during zebrafish embryo development" (see preprint doi:10.11012024.01.30.575011).


## setup

To rerun the analysis, clone this git repository, install dependencies, copy input files to `raw` folder, set up sample meta data file and edit config file as needed.

### install dependencies with conda

After installing conda (see: https://docs.anaconda.com/free/miniconda/miniconda-install), the required dependencies e.g. can be easily installed into a conda enviornment using the conda create:

```
$ conda create --name codonusage \
  bioconda::snakemake \
  anaconda::pandas \
  conda-forge::scikit-learn \
  conda-forge::numpy \
  conda-forge::matplotlib-base \
  conda-forge::seaborn-base \
  conda-forge::biopython
$ conda activate codonusage
```

Note: On windows snakemakecan not be installed though bioconda, other means are available. Please see https://snakemake.readthedocs.io/en/stable/getting_started/installation.html.

### set up input files

Input files should be provided within the `raw` folder:

- `gene_seqeunces`: fasta file with all zebrafish coding sequences as provided by ensembl (https://ftp.ensembl.org/pub/release-[releasID]/fasta/danio_rerio/cds/; adjust releaseID to fit your gene_counts file and appris annotation)

- `gene_counts`: file with gene counts per sample. It contains one row per gene; one column stating the `GeneType`, eg "protein coding" or "rRNA"; one column stating `GeneSymbol`, one column with the ensembl gene id (`ENSG`) and one column per sample respresenting the read count per gene.

- `apris_data`: isoform annotation from appris (https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/danio_rerio/GRCz11/appris_data.appris.txt [as 16th of June 2024: ensembl104]

### set up sample meta data

Sample meta data should be provided in a tsv file. By default it is located in the `config` folder and named `samples.tsv`.
The tsv file has to include 3 columns with corresponding headers:

- `sample`: containing the sample id;

- `timepoint`: the 1-based index of the timepoint within the analysed time series at witch the samples was draw;

- `timepiont_name`: text description of the time point.

An example is shown below:

```
sample	timepoint	timepoint_name
EV06001	1	ovary
EV06008	1	ovary
EV09001	1	ovary
EV06002	2	eggs
EV06009	2	eggs
EV09002	2	eggs
EV06003	3	4-cell, 1 hpf
EV06010	3	4-cell, 1 hpf
...
```

### set up config

The input files and some parameters are specified in the `config.yaml` file found in `config` folder.

- `gene_isoforms_fasta`: 'Danio_rerio.GRCz11.cds.all.fa' -- name of fasta file with all zebrafish coding sequences

- `appris_annotation`: 'appris_data.appris.txt' -- name of appris annotation file

- `samples_tsv`: 'config/samples.tsv' -- path to file with sample meta data as described in "set up sample meta data" section

- `raw_sequencing_counts`: 'raw_gene_counts_GeneSymbols_GeneTypes.tsv' -- name of gene counts file

- `isoform_selection`: 'priority_principal_triple_longest' -- method for selecting representative isoform of gene

- `RPM_cutoff`: .000010 -- minimum coverage of a gene to be included in analyses

```
RPM        fraction
0000000    .000000
0000001    .000001
0000010    .000010
0000100    .000100
0001000    .001000
1000000    1
```

## run

```
snakemake 'results/all_done.txt' -n 1

```

