# GeneNetwork Pipeline

## Introduction
The GeneNetwork pipeline is a command line python3 tool to perform the decomposition and gene pathway score calculation. The main goal of the analysis is to calculate gene z-scores for every gene for every pathway in the input files which gives the contribution of a gene to the pathway. Where the input pathway file contains only a Boolean value if the gene is linked to the pathway, this analysis calculates a contribution z-score value for each gene to every pathway by using gene expression data.

The entire method and clinical application of this pipeline is described in our article “KidneyNetwork: Using kidney derived gene expression data to predict and prioritize known and novel genes involved in kidney disease” (under review). Below you find a summary of the method and a description how the pipeline can be used.

## Pipeline description
The pipelines contain two steps: 1) the decomposition and 2) the gene-phenotype score calculation. The two steps are separated in two different python scripts.

## Decomposition
The first step of the pipeline is the calculation of the eigenvectors from the gene expression dataset. Three different decomposition methods are implemented, PCA, FastICA and a stabilized version of FastICA. The decomposition method PCA was implemented by using Singular value decomposition class of Sklearn package and the svd package of SciPy. The different algorithms van be selected by using the “pca_type” parameter. See parameter description for more information.

The input of the script in a tab separated gene expression file where the rows represented the samples and the columns the genes. The output contains one file with the gene/component’s matrix (eigenvectors) and one file with the sample/component’s matrix (projected data). If the “input_transpose” option was set, both output files are swapped to maintain the structure and layout of both files. Beside the output matrixes a log file is created.

The output is by default exported to both txt and pickle format. The .txt file is gzip compressed. If some files are not necessary, they can be disabled by using the output parameters.

## Gene-phenotype score calculation
The second step is the calculation of the gene-pathway scores based on the eigenvectors obtained from the gene expression data and the give(n) pathway databases. The given eigenvectors were used to fit a logistic model to identify an intercept and the coefficients which were used to calculate a gene log_odds score. This score is translated to z-scores by using a permutation strategy. The gene-components matrix is permutated to create a null-distribution. The fitted intercept and coeffects were used to calculate the gne log_odds scores corresponded to the null distribution and the mean and standard deviation (std) of these values were calculated. The mean and std are used to translate the gene log_odds score to z_scores.

The final gene z-scores were used to calculation the performance of the model by calculating the Area under the ROC-curve (AUC value) and the corresponded p-value by using the Mann-Whitney rank test. The p_values are also Bonferroni corrected.

The input files of the second step of the pipeline consist of the eigenvector matrix from created by the decomposition step. The input can be both the .txt file or the pickle file. The second input file contains the pathway matrix file which indications which gene is linked to which pathway. The file must contain a tab separated matrix where the rows correspond to the pathways and the columns to the genes. The values must be Booleans indicating is a gene is linked (value is 1) or not linked (value is 0) to a pathway.

The model is by default fitted on the complete gene list, but a gene sub selection can be given to fit the model only on a subset of genes. The subset can be given as path to a file with the gene names, where every row corresponds to one gene.

The output of the pipeline contains one file with the gene z-scores for every pathway (gene_pathway_scores) where the rows correspond to every pathway and the columns to every gene. The second output file contains the prediction AUC values and the corresponded p value and Bonferroni corrected p value. The output is by default exported to both txt and pickle format. The .txt file is gzip compressed. If some files are not necessary, they can be disabled by using the output parameters.

The permutated version of the input matrix will be saved as pandas pickle file on the same location of the input file with the suffix “_permutation_matrix”. If this matrix is already present, the script will read this file instead of creating a new permutated matrix (except when the option –force is set)

The process of calculation of z-scores for one pathway can be take a serious amount of time, especially for pathways with a lot of genes. To reduce the overall processing time, the script use by default multiple CPU cores of the machine. The maximal number of available cores are used, but this can be reduced by setting a maximal number of cores through a parameter.  The script has also build-in support for multi node processing, for cluster usage. This can be easily setup with a job manager like Slurm. An example job script is included.

## Input data
Both scripts use input data frames to perform the analysis. Every input data frame can be given as Pandas pickled data frame or as a tab separated txt file. The txt input file is automatically decompressed if the file paths ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’. An cached pandas pickle file are created on the same location as the input txt file to speed up the reading process of the input file in further analysis. These cached files automatically created on the same location and same name as the input file with the suffix “_cashed.pickle”.

## Output data
Both scripts create different output files as described in the above. The files are created in the given output directory, which is also created if the directory does not exit. Output matrixes are by default created as pandas pickled dataframe and tab separated text file. The txt matrix files are by default gz compressed. The compressing can be disabled (parameter: --output_disable_gzip), even as the creation of the txt file or the pickle file (option: --output_disable_txt and --output_disable_pickle)

# Prerequisites  

This program is developed in Python 3.7, performance on other versions is not guaranteed.

The program requires the following packages to be installed: 

* pandas ([v0.25.1](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))
* numpy ([v1.17.2](https://github.com/numpy/numpy/releases); [BSD License](https://www.numpy.org/license.html))  
* sklearn ([v0.22.2](https://scikit-learn.org/0.22/); [BSD 3-Clause License](https://github.com/scikit-learn/scikit-learn/blob/master/COPYING))    
* matplotlib ([v3.1.0](https://github.com/matplotlib/matplotlib/releases); [PSF License](https://matplotlib.org/3.1.0/users/license.html))
* seaborn ([v0.9.0](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE)) 

See 'Installing' on how to install these packages.

## Installing  

Installing and running the program can be done by executing the following steps:  

**Step 1: acquire the source files**
Download or clone the GeneNetworkAnalysis directory from the [systemsgenetics](https://github.com/molgenis/systemsgenetics) repository.

**Step 2: installing the required packages**  
1. Install the required packages. The python packager manager (pip) and the requirements file can be used to install all the nessesary packages:  
    ```  
    pip install -r requirements.txt
    ```  
2. Wait for the command to finish

# Scripts

De repository contains different scripts to create and analyse gene expression networks.
Most of the scripts contain an argument parser and help function.

## Perform decomposition
The decomposition step is implemented in the “gene_network_decomposition.py” script. The help and parameter instructions can be found by using the command:

    ```  
    python3 gene_network_decomposition.py --help
    ```  

### Parameters
| Short     | Long                                             | Description                                                                                                                                                                                                                                                                                                                                                 |
|-----------|--------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| &#x2011;i | &#x2011;&#x2011;input                            | The path of the input file. The input file contains the gene expression values where the rows correspond to the samples and the columns to the genes. The input file can be given as pandas pickled dataframe or as tab separated txt file. The txt input file is automatically decompressed if the file paths ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’.  |
| &#x2011;o | &#x2011;&#x2011;output                           | Path to the output directory. The output files will be created in this directory. The directory will automatically be created if not exists.                                                                                                                                                                                                                |
| &#x2011;a | &#x2011;&#x2011;analysis                         | Set the analysis type which will be used to perform the decomposition. Can be set to ‘PCA’, ‘FASTICA’, ‘FASTICA_STABLE’ or ‘FASTICA_STABLE_FROM_TEMP’ PCA: Use principal component analysis as decomposition method. <br><br> PCA is implemented by using SVD algorithms Different algorithms are implemented and can be selected by using the ‘pca_type’ parameter. <br><br> FASTICA: Use the fastICA implementation to perform an independent component analysis as decomposition method. The whiting step is done by the implemented SVD algorithms and can be selected by using the ‘pca_type’ parameter. <br><br> FASTICA_STABLE: Use multiple different fastICA runs and calculate the average gene eigenvector loadings of the independent components to get sable independent components. The FastICA algorithms use a random initialization which can result in slightly different outcomes per individual run. This analysis can take a long time. <br><br> FASTICA_STABLE_FROM_TEMP: Method to finish up an uncomplete “FASTICA_STABLE” run. The executed runs will be combined to the final independent components. The method will not execute more FastCIA runs even if the number of set iterations is not reached (the option: fastica_stable_iterations) |
| &#x2011;t | &#x2011;&#x2011;test                             | Perform a test run with the first 150 samples and 100 genes                                                                                                                                                                                                                                                                                                 |
| &#x2011;c | &#x2011;&#x2011;components                       | The number of components which will be returned in the export files. If the number is not set, the maximal number of components will be returned                                                                                                                                                                                                            |
|           | &#x2011;&#x2011;fastica_stable_iterations        | Number of FastICA runs used to calculate the final FastICA stable results. Can only be set if analysis is set to FastICA-stable. Default is 10                                                                                                                                                                                                              |
|           | &#x2011;&#x2011;fastica_stable_safe_intermediate | Safe the intermediate FastICA runs in the fastICA-stable run. Can only be set if analysis is set to FastICA-stable.                                                                                                                                                                                                                                         |
| &#x2011;n | &#x2011;&#x2011;nrows                            | Number of rows in the input dataset to process.                                                                                                                                                                                                                                                                                                             |
|           | &#x2011;&#x2011;input_transpose                  | Transpose the input matrix before the analysis. The columns will be processed as rows and the rows as columns. The output matrixes are not changed.                                                                                                                                                                                                         |
|           | &#x2011;&#x2011;force                            | Do not use the cached versions of the input data but the orginal input files.                                                                                                                                                                                                                                                                               |
|           | &#x2011;&#x2011;log2                             | Perform log2 transformation of the input data before the analysis                                                                                                                                                                                                                                                                                           |
|           | &#x2011;&#x2011;center_scale                     | Perform center scale transformation of the input data before the analysis                                                                                                                                                                                                                                                                                           |
|           | &#x2011;&#x2011;output_disable_gzip              | Do not gzip the data output text files                                                                                                                                                                                                                                                                                                                      |
|           | &#x2011;&#x2011;output_disable_txt               | Do not export the output as text files.                                                                                                                                                                                                                                                                                                                     |
|           | &#x2011;&#x2011;output_disable_pickle            | Do not export the output as pickle Pandas files.                                                                                                                                                                                                                                                                                                            |


## Perform gene pathway score prediction
The decomposition step is implemented in the “gene_pathway_score_prediction.py” script. The help and parameter instructions can be found by using the command:

    ```  
    python3 gene_pathway_score_prediction.py --help
    ```  
### Parameters

| Short     | Long                                  | Description                                                                                                                                                                                                                                                                                                |
|-----------|---------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| &#x2011;i | &#x2011;&#x2011;input                 | The path to the input components / genes matrix. The input file can be given as pandas pickled dataframe or as tab separated txt file. The txt input file is automatically decompressed if the file path ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’.                                                       |
| &#x2011;o | &#x2011;&#x2011;output                | Path to the output directory. The output files will be created in this directory. The directory will automatically be created if not exists.                                                                                                                                                               |
| &#x2011;p | &#x2011;&#x2011;pathway_matrix        | The path to the pathway matrix file(s). The parameter must contain at least one path to a pathway matrix but can also contain multiple matrixes by separating the paths by commas. <br> The input file(s) can be given as pandas pickled dataframe(s) or as tab separated txt file(s). The txt input file is automatically decompressed if the file path ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’.                                                                                               |
| &#x2011;g | &#x2011;&#x2011;genes                 | The path to the background gene file(s). The genes mentioned in the file(s) are used to fit the model instead of all the genes present in the input matrix. If multiple pathway matrixes were given, than the background files must be given in the same order as the pathway matrix. <br> The parameter must contain at least one path to a pathway matrix but can also contain multiple matrixes by separating the paths by commas. <br> The input file(s) can be given as pandas pickled dataframe(s) or as tab separated txt file(s). The txt input file is automatically decompressed if the file path ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’.                                                                                               |
| &#x2011;a | &#x2011;&#x2011;analysis              | Set the analysis which will be used to fit the model and calculate the gene pathway scores. Available analysis types: “regression” and “t_test”. <br><br> regression: The updated method which used the logistic regression model to fit an model which is used to calculate the gene z_scores by using permutation of the input matrix. <br><br> t_test: The original gene network method which use an t-test apparoach combined with an linear model to calculate the gene p-values which are translated to the gene z-scores                                                                                                                              |
| &#x2011;m | &#x2011;&#x2011;minimal_genes         | The minimal number of genes which must be present in the pathway to be analyzed. Pathways with a lower number of genes are skipped.                                                                                                                                                                        |
|           | &#x2011;&#x2011;start                 | The first index of the pathway matrix to start the analysis. Can be used to subset the input pathway matrix (together with --end).                                                                                                                                                                         |
|           | &#x2011;&#x2011;end                   | The last index of the pathway matrix to end the analysis. Can be used to subset the input pathway matrix (together with --end).                                                                                                                                                                            |
|           | &#x2011;&#x2011;node                  | The id of the single node (zero first). This option can be used to execute the script in an multi node environment (cluster). Every node will execute an equal number of pathways. Option must be used together with the parameter “num_nodes”. <br><br> For example, if the cluster has 10 nodes and the pathway file contains 100 pathways, every node will process 10 pathways. The first node (--node 1) will process the first 10 pathways, the second node (--node 2) will process pathway 11 until 20, the third node will process pathway 21 until 30, etc. |
|           | &#x2011;&#x2011;num_nodes             | Total number of nodes used in the cluster (multi node) environment. This option can be used to execute the script in an multi node environment (cluster). Every node will execute an equal number of pathways.  |
| &#x2011;c | &#x2011;&#x2011;cores                 | Number of cores which will be used to process the pathways on different processor cores. Default is the number of available.                                                                                                                                                                               |
|           | &#x2011;&#x2011;force                 | Do not use the cached versions of the input data but the original input files.                                                                                                                                                                                                                             |
|           | &#x2011;&#x2011;output_disable_gzip   | Do not gzip the data output text files                                                                                                                                                                                                                                                                     |
|           | &#x2011;&#x2011;output_disable_txt    | Do not export the output as text files.                                                                                                                                                                                                                                                                    |
|           | &#x2011;&#x2011;output_disable_pickle | Do not export the output as pickle Pandas files.                                                                                                                                                                                                                                                           |

## Author  
Henry Wiersma *(1)*  
  
1. Genetics Department of the University Medical Center Groningen (UMCG), The Netherlands.
