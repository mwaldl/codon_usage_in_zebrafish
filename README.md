# analysing codon usage in zebrafish

pipeline used to analyse codon usage during zebrafish embryo development as presented in "tRNA expression and modification landscapes, and their dynamics during zebrafish embryo development" (see preprint doi:10.11012024.01.30.575011).


## setup

To rerun the analysis, clone this git repository, install dependencies, copy input files to `raw` folder, set up sample meta data file and edit config file as needed.

### clone workflow
Download this workflow, e.g. with git clone:

`git clone git@github.com:mwaldl/codon_usage_in_zebrafish.git`

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

Input files should be provided within the following sub-folders of the  `raw` folder:

- `gene_seqeunces`: a fasta file with all zebrafish coding sequences as provided by ensembl (https://ftp.ensembl.org/pub/release-[releasID]/fasta/danio_rerio/cds/; adjust releaseID to fit your gene_counts file and appris annotation)

- `gene_counts`: a file with gene counts per sample. It contains one row per gene; one column stating the `GeneType`, eg "protein coding" or "rRNA"; one column stating `GeneSymbol`, one column with the ensembl gene id (`ENSG`) and one column per sample respresenting the read count per gene.

- `apris_data`: a file with isoform annotation from appris (https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/danio_rerio/GRCz11/appris_data.appris.txt [as 16th of June 2024: ensembl104]

### set up config

The input files and some parameters are specified in the `config.yaml` file found in `config` folder. 

- `gene_isoforms_fasta`: 'Danio_rerio.GRCz11.cds.all.fa' -- name of fasta file with all zebrafish coding sequences located with in the raw/gene_seqeunces folder

- `appris_annotation`: 'appris_data.appris.txt' -- name of appris annotation file located with in the raw/apris_data folder

- `samples_tsv`: 'config/samples.tsv' -- path to file with sample meta data as described in "set up sample meta data" section

- `raw_sequencing_counts`: 'raw_gene_counts_GeneSymbols_GeneTypes.tsv' -- name of gene counts file located with in the raw/gene_counts folder

- `isoform_selection`: 'priority_principal_triple_longest' -- method for selecting representative isoform of gene

- `RPM_cutoff`: 0 -- minimum coverage in RPM of a gene to be included in analyses

### set up sample meta data

Sample meta data should be provided in a tsv file. By default it is located in the `config` folder and named `samples.tsv`.
The tsv file has to include 3 columns with corresponding headers:

- `sample`: containing the sample id;

- `timepoint`: the 1-based index of the time point within the analysed time series at witch the samples was draw;

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

## run the pipeline

1. Activate the newly set up conda environment if it is not yet activated:
   ```
   conda activate codonusage

   ```

2. Navigate to the workflow folder:
   ```
   cd codon_usage_in_zebrafish

   ```

3. Run the pipeline:
   ```
   snakemake 'results/all_done.txt' -c 2

   ```
   replace '2' with wathever number of cores you want to use

## output

### codon usage per time point
Codon usage per time point and per sample is provided as tsv file and as heatmap in the results/codon_usage_per-timepoint folder. The file name consists of four parts: a descriptor for the method used to select the principle isoform, wether the data is provided per sample or as mean over three samples per time point, wether the rows ar sorted and labled by codon or by the name of the decoded amino accid and if the file contains the data table or the heatmap. `<isoform selection method>_isoform_<time/sample>_vs_<codon/aminoacid_codon>.<pdf/tsv>`

### comaprison of isoform selection methods
The folder results/method_comparison contains two heatmaps that show the mean codon fraction over all samples coputed with different isoform selection methods. The file 'method_cu_fraction.pdf' shows the fraction that each codon contributes. The file 'method_cu_count.pdf' shows the rsult of multiplying the codon count in each gene with the abundance of each gene. The first plots shows no major differnces between methods, while the second methods shows variations in the total number of codons. The correspodning data files can be found in the resources/method_comparisonfolder.

### pca

PCA was performed for gene expression levels as well as codon fractions. 

The PCA on gene expression levels only considers genomic protein coding genes (see tsv file in resources/normalized_counts folder). It is performed for per sample gene expression data and on the mean gene expression level per time point. Accordingly, the results are found in results/pca/gene_expression/sample and results/pca/gene_expression/timepoint. A two 2D plot of the first two components is shown in pca.pdf; the explained variance per component is plotted in explained_variance.pdf and the genes that represent at least 10% to a component are shown in explained_variance.pdf_feature_contribution_bar.pdf, where the first plot represents the first component, second plot the second component etc. 

The PCA on codon usage is located in the results/pca/codon_usage folder, seperated by subfolders based on the isoform selction method. The pca.pdf, explained_variance.pdf and explained_variance.pdf_feature_contribution_bar.pdf files are produced analogously to the gene expression PCA plots. Contributions to PCs are of course  codons fractions rather than gene expression fractions.


