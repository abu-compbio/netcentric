# A Network-centric Framework for the Evaluation of Mutual Exclusivity Tests on Cancer Drivers

This is the original repository for the Network-centric Framework for the Evaluation of Mutual Exclusivity Tests project codes. The project involves the evaluation of mutual exclusivity methods: Discover, Discover_strat, Fisher's Exact Test, MEGSA, MEMO, and WExT. The results from these methods include pairwise mutual exclusivity p-values. Based on them, we apply our network-centric epistatic evaluation.


## Getting Started
### Prerequisites
Python: 3.8-3.9

### Getting Started

Using Github clone

```bash
git clone  https://github.com/abu-compbio/netcentric
cd NetCentric
pip install -r requirements.txt
```
You can also run this project locally by following these steps:
1. Download the repo
2. Unzip NetCentric-main
3. Open cmd/terminal and cd into the project
4. Execute python -m pip install -r requirements.txt

Some of the datasets are given in the .zip file format. In order to unzip, you should run script_unzip_data.py located under the NetCentric.

```bash
python script_unzip_data.py
```

## Input

### 1. The PPI network edge and index files. 

The file is located at data/intact_nodupl_index_file.txt where the first column contains gene identifiers and the second column corresponds to gene symbols.

```bash
0	""CHEBI
1	100147744
2	1B
3	1EFV
...
``` 
The file is located at data/intact_nodupl_edge_file.txt where the first two columns contain the gene identifiers and the third column denotes the confidence level of the interaction. 

```bash
7589	13441	0.99
9123	10446	0.98
4248	1740	0.98
3776	5279	0.98
...
``` 

### 2. Input mutation data

The mutation data includes pairwise mutual exclusivity p-values given for each method (discover, discover_strat, fishers, megsa, memo and wext).
The files with the name of mutations_all_genes include all genes and intact_filtered include only ones in intact network. 

The file is located at data/{method}_mutation_filtered_ep_data/{cancerType}_{method}_result_mutations_all_genes_{threshold}.txt
```bash
	gene1	gene2	pvalue
0	A2M	A2ML1	0.6584654889330113
1	A2M	ABCA1	0.5332913581418495
2	A2M	ABCA10	0.8971732886956303
...
``` 
The file is located at data/{method}_mutation_filtered_ep_data/{cancerType}_{method}_pairs_intact_filtered_subset{threshold}.txt

```bash
	gene1	gene2	pvalue	oddsratio
0	TCF7L2	CTNNB1	0.9015805073650888	1.6786858974358974
1	SMAD4	SMAD3	0.839665475908354	1.4567901234567902
2	EP300	TP53	0.0742406168447767	0.5221052631578947
...
``` 

The file is located at data/binary_matrices_all_genes_ep_mutation_filtered/ directory. Each row is a TCGA patient id and each column is a gene. The matrix contains 1 if the gene is mutated in the corresponding patient. Here, we only provide the mutation matrix  for COADREAD.

```bash
	patients	A1BG	A1CF	A2M ...
TCGA-3L-AA1B-01A	0	0	0
TCGA-4N-A93T-01A	0	0	0
TCGA-4T-AA8H-01A	0	0	0
...
``` 

MLA: The file contains the corresponding MLA.

The file is located at data/MLA_ep_mutation_filtered_all_genes

```bash
A1BG	4.261253658699028
A1CF	5.095042391780406
A2M	5.539871662596874
...
``` 

### 3. Reference cancer genes

The file is located at data/known_cancer_genes directory.

 1. CGC genes:
We download all the genes from Cancer Gene Census from COSMIC database.

 2. CGC_SNV genes:
We try using a subset of CGC genes to include only those which have SNV type of mutations in cancer (378 out of 723 genes). To this end, we filter out the genes where the mutation type column consists of only A (amplification), D (large deletion) or T (translocation). 

 3. IntoGen genes:
We download Unfiltered driver results 05.tsv file (2020-02- 02 release) from https://www.intogen.org and include the genes where FILTER column is PASS, which results in 503 genes.  


### 4. TSN

The file contains gtex edges for the corresponding tissue type with a given threshold (0.0 and 0.5)

The file is located at data/gtex_tsn_fractions_intact_filtered_applied_threshold

```bash
MDM2	TP53	1.0
PAK1	RAC1	1.0
FADD	CASP8	0.9987163029525032
...
``` 

### 5. Results

The Mutual Exclusivity results will be available in the folder ME_results.
The TSN results will be available in the folder tsn_results.
The MLA results will be available in the folder MLA_results.
These folders will appear underthe main directory, when the results are ready.


## Runs

The codes regarding various analyses given in the main article.

### **ME Evaluations Based on Defined Metrics** 

The main source code for the evaluation ME Tests. In the main article it was discussed under the section "ME Evaluations Based on Defined Metrics".
This analysis code were also used in the section "Robustness Analysis of Evaluations Based on Defined Metrics". robustness_iterations value is given as parameter i: number of iteration in the code. 

As output, you get tables with all analysis results in NetCentric/ME_results
To generate the algorithm for the given input, the following script should be run 
(c: cancer type, t: threshold, i: number of iteration, m: methods, p: p_value threshold, -ni: network index file, -e:network edge file, -r: Reference cancer genes)

```bash
cd src
evaluations_on_metrics.py -c COADREAD -t 20 -i 100 -m discover discover_strat fishers megsa memo wext -p 0.05 -ni intact_nodupl_index_file.txt -e intact_nodupl_edge_file.txt -r Census_allFri_Apr_26_12_49_57_2019.tsv
``` 

### **ME Evaluations Based on Corrections via MLA**

Scatterplots of percentage significance of mutual exclusivity runs vs mutation load association (MLA). In the main article it was discussed under the section "ME Evaluations Based on Corrections via MLA". As output, you get results in NetCentric/MLA_results/percent_sig_figures

```bash
cd src
evaluations_via_mla.py -c COADREAD -t 20 -m discover discover_strat fishers megsa memo wext
```

Scatterplots of percentage significance of mutual exclusivity runs vs mutation load association (MLA) when only CGC genes that have > 1 neighbors are included.
As output, you get results in NetCentric/MLA_results/perc_sig_figures_for_multiple_neighbors

```bash
cd src
evaluations_via_mla_neighbors.py -c COADREAD -t 20 -m discover discover_strat fishers megsa memo wext
```

### **ME Evaluations Based on Corrections via TSN**

In the main article it was discussed under the section "ME Evaluations Based on Corrections via TSN".
(c: cancer type, t: threshold, m: methods, ti: tissue, th: tsn threshold)
As output, you get tables with all analysis results in NetCentric/tsn_results.

```bash
cd src
evaluations_via_tsn.py -c COADREAD -t 20 -m discover discover_strat fishers megsa memo wext -ti Colon -th 0.0
```

ROC curves for comparing the mutual exclusivities of tissue-specific and non-tissue specific CGC-CGC gene pairs and non-CGC-non-CGC gene pairs on.
As output, you get results in NetCentric/tsn_results/figure_tsn_AUROC
(c: cancer type, t: threshold, m: methods, th: tsn threshold, p: percentage )

```bash
cd src
me_on_tsn_ntsn_roc_curve.py -c COADREAD -t 20 -m discover discover_strat fishers megsa memo wext -th 0.0 -p 0.25

```

