# Following are the instructions to run the mutual exclusivity methods with to get the desired results:

## 1. DISCOVER (Canisius et. al. 2016)
#### Github repository: https://github.com/NKI-CCB/DISCOVER
We used the python version of discover for Netcentric. At first create a conda environment by using the following command:

```
conda create -n discover -c http://ccb.nki.nl/software/discover/repos/conda discover
```

After creating and activating the environment, follow the instructions in the documentation: https://ccb.nki.nl/software/discover/doc/python/discover-intro.html

Make sure your mutation matrix follows the input format shown in the link above. Then simply load the data and apply the pairwise test command on it

```python
mut = discover.datasets.load_dataset("data")
events = discover.DiscoverMatrix(mut)
subset = mut.sum(1) > threshold
result_mutex = discover.pairwise_discover_test(events[subset])
```

In the above script, data is the data you provide and threshold is the mimimum mutations allowed. The result is a table of p-values and q-values that looks like the following:

```
     gene1   gene2        pvalue        qvalue
0   ARID1A    TP53  1.707299e-04  3.720345e-03
1     CDH1    FAT3  2.367740e-04  4.229149e-03
```

-----------------------
 
## 2. DISCOVER_strat
This method was only applied in the case of BRCA and COADREAD cancer types where the samples were clearly stratified among well defined groups. 

In this case, first generate a stratified list that contains the groups that the samples belong to (same order). Then simply change events to the following:

```
events = discover.DiscoverMatrix(df2, strata=stratified_list)
```

The rest should be run as before and the results will also be similar. For more information follow the tutorial: https://ccb.nki.nl/software/discover/doc/python/discover-intro.html

-----------------------

## Fisher's Exact Test
The Fisher's exact test was performed by simply using scipy.stats.fisher_exact command from the scipy library. In order to perform the test between two genes A and B, consider the samples where A and B are mutated. Then calculate the following events:

I. both_A_and_B: samples where both gene A and gene B mutate

II. B_but_not_A: samples where gene B mutates but not A

III. A_but_not_B: samples where gene A mutates but not B

IV. no_A_B: all samples - union(samples_A, samples_B)

Then the pairwise Fisher's exact test is performed using the following command:

```python
oddsratio, pvalue = fisher_exact([[both_A_and_B, B_but_not_A],[A_but_not_B,no_A_B]], alternative='less')
```

-----------------------

## 4. MEMO (Ciriello et. al. 2011)
#### Package available: http://ciriellolab.org/dataandtools.html
This method was developed in Java by the authors. However it proved difficult to implement and adjust with our data. We rewrote the code in python and ran it on our data.

The "MEMO_bipartite_from_binary_graph_alternate" notebook contains the functions and step by step instructions used for running the pyhton version of MEMO.

-----------------------

## 4. MEGSA (Hua et. al. 2016)
#### R Package available: https://dceg.cancer.gov/tools/analysis/megsa
MEGSA is a statistical test run on R. From the main script of MEGSA, the funestimate function was used to get the log likelyhood between two genes. The log likelyhood value was then used to perform a chi square test that returned the desired pairwise p-values. Follow the adjusted script, "run_megsa.R" for the commands to run MEGSA.

It is to be noted that, MEGSA authors discouraged running MEGSA on a dataset of more than 100 genes. Therefore, MEGSA runs were limited to the $t=20$ setting.

-----------------------
  
## 4. WExT (Leiserson et. al. 2016)
#### Github repository: https://github.com/raphael-group/wext
For WExT follow the instructions given in their github repository and run "find_exclusive_sets.py" script with WRE Saddlepoint approximation. The gene set size, k, was limited to 2 for pairwise results.wext. Combining all the parameters, the following commands were used:

1. Define the parameters:  
```
num_permutations=10000
num_cores=30
C=COADREAD
T=5
```

2. Compute mutation probabilities:  
```
python2 ../../compute_mutation_probabilities.py \
	-mf adjacency_files/${C}_adjacency_json_t${T}.json \
	-np $num_permutations \
	-nc $num_cores \
	-wf weights_files/${C}_weights_t${T}.npy \
	-s  12345 \
	-v  1
```

3. Find sets using mutual exclusivity test statistic:  
```
python ../../find_exclusive_sets.py \
	-mf adjacency_files/${C}_adjacency_json_t${T}.json \
	-ks  2 \
	-c  $num_cores \
	-f  0 \
	-s	Enumerate \
	-o  exclusivity_results/${C}_wext_raw_output_t${T} \
	-v  4 \
	WRE \
	-m Saddlepoint \
	-wf weights_files/${C}_weights_t${T}.npy \
```

Note that the adjacency file contains an adjacency matrix that follows the input format for WExT.