HotNet
=======================

**We have introduced a new method, Hierarchical HotNet, that improves on HotNet.  As a result, HotNet is no longer actively updated.  Please see the Hierarchical HotNet [manuscript](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236) and [GitHub repository](https://github.com/raphael-group/hierarchical-hotnet) for details.** 

___

HotNet is an algorithm for finding significantly altered subnetworks in a large gene interaction
network. While originally developed for use with cancer mutation data, the current release of
HotNet also supports the application in which scores can be assigned to genes in the network
(previously called generalizedHotNet).


Requirements
------------------------

* Linux/Unix
* [Python 2.7](http://python.org/)
* [NumPy 1.6.2](http://www.numpy.org/)
* [SciPy 0.10.1](http://www.scipy.org/)
* [NetworkX 1.7](http://networkx.github.io/)

HotNet will likely work with additional versions of Python, NetworkX, NumPy, and SciPy, but
alternative configurations have not been tested.

Support
------------------------
For support using HotNet, please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users).

Simple runs
------------------------
To get started running HotNet quickly and easily, use the `simpleRun.py` Python script.  You must
provide the following parameters:

        =====================================================================================
        | PARAMETER NAME         | DESCRIPTION                                              |
        =====================================================================================
        |-mf/--infmat_file       |Path to .mat file containing influence matrix downloaded  |
        |                        |from http://compbio.cs.brown.edu/projects/hotnet/         |
        -------------------------------------------------------------------------------------
        |-if/--infmat_index_file |Path to gene-index mapping file downloaded from           |
        |                        |http://compbio.cs.brown.edu/projects/hotnet/              |
        -------------------------------------------------------------------------------------
        |-hf/--heat_file         |Path to a tab-separated file containing a gene name in the|
        |                        |first column and the heat score for that gene in the      |
        |                        |second column of each line.                               |
        -------------------------------------------------------------------------------------
        
Running with only the parameters specified above will create a 'hotnet_output' directory in your
current working directory that contains 5 subdirectories each prefixed with `delta_`. Each of these
subdirectories contains results files for a different value of the delta parameter used by the
HotNet algorithm. The output files are:

* `components.txt`: Lists subnetworks identified as significantly altered, one per line. Genes
  within each subnetwork are separated by tabs.
* `significance.txt`: For k from 2 - 10, lists the number of subnetworks of size >= k found in the
  real data, the expected number of subnetworks of size >= k based on permuted data, and the p-value
  for the observed number of subnetworks.
* `results.json`: Contains all of the above information plus the parameters used for the run in
  JSON format to faciliate further automated processing
* `heat.json`: Input heat scores in JSON format to facilitate further automated processing.

The `simpleRun.py` script can also be used to create a web visualization of the output subnetworks.
To do so, include the `--edge_file` parameter:

        ========================================================================================================
        | PARAMETER NAME         | DEFAULT          | DESCRIPTION                                              |
        ========================================================================================================
        |-ef/--edge_file         | None             |Path to TSV file listing edges of the interaction network,|
        |                        |                  |where each row contains the indices of two genes that are |
        |                        |                  |connected in the network. This is used to create          |
        |                        |                  |subnetwork visualizations; if not provided, visualizations|
        |                        |                  |will not be made.                                         |
        --------------------------------------------------------------------------------------------------------
        |-nn/--network_name      | Network          |Display name for the interaction network.                 |
        --------------------------------------------------------------------------------------------------------

This will result in a a `viz` subdirectory of the output directory. To view the visualizations,
navigate to the `viz` directory and run `python -m SimpleHTTPServer`, then visit `http://localhost:8000`
in a browser.

To see an example, first make sure you have downloaded the influence matrices from
[http://compbio.cs.brown.edu/projects/hotnet/](http://compbio.cs.brown.edu/projects/hotnet/)
and saved them in the `influence_matrices` directory, then run:

    python simpleRun.py @example/configs/simple.config

When using `simpleRun.py`, you may also optionally provide any or all of the parameters listed
below. If one of these parameters is not provided, it will be set to the default value shown below.

        ========================================================================================================
        | PARAMETER NAME         | DEFAULT          | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size      | 3                |Minimum size connected components that should be returned.|
        --------------------------------------------------------------------------------------------------------
        |-ms/--min_heat_score    | See description  |Minimum heat score for a gene to be eligible for inclusion|
        |                        |                  |in a returned connected component. By default, all genes  |
        |                        |                  |with positive heat scores will be included. (To include   |
        |                        |                  |genes with score zero, set min_heat_score to 0).          |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | 100              |Number of permutations that should be used for parameter  |
        |                        |                  |selection and statistical significance testing            |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | hotnet_output    |Output directory.                                         |
        --------------------------------------------------------------------------------------------------------


Advanced use
------------------------
The `simpleRun.py` script described above runs the entire HotNet pipeline in one command. For more
advanced use cases, you can also perform each step individually. In particular, you may wish to
follow the steps below when using mutation data.

The steps of the Hotnet algorithm and the code provided for each step are described below.

1. ###Influence matrix creation###

    This step creates a matrix that defines an "influence score" for each gene pair in the network
    based on known gene interactions and a heat diffusion process. Code is not provided for
    influence matrix creation. Instead, please use one of the pre-computed influence matrices
    provided.


2. ###Heat score generation###

    This step creates a JSON file containing heat scores on each gene required in subsequent steps.
    Heat scores can either be specified directly or calculated 

    The Python script `generateHeat.py` can be used to create the JSON heat file. The required and
    optional parameters to the script are described below.

    If heat scores are specified directly, the first argument to `generateHeat.py` should be
    `scores`, e.g.:

            python generateHeat.py scores <additional_parameters>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |--heat_file             | REQUIRED         |Path to a tab-separated file containing a gene name in the|
        |                        |                  |first column and the heat score for that gene in the      |
        |                        |                  |second column of each line.                               |
        --------------------------------------------------------------------------------------------------------
        |--gene_filter_file      | None             |Path to file listing genes whose heat scores should be    |
        |                        |                  |preserved. If present, heat scores for all genes not      |
        |                        |                  |listed in the file will be discarded.                     |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------

    If heat scores are to be calculated from mutation data, the first argument to `generateHeat.py`
    should be `mutation`, e.g.:

            python generateHeat.py mutation <additional_parameters>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |--snv_file              | REQUIRED         |Path to a tab-separated file containing single nucleotide |
        |                        |                  |variants(SNVs) where the first column of each line is a   |
        |                        |                  |sample ID and subsequent columns contain the names of     |
        |                        |                  |genes with SNVs in that sample. Lines starting with '#'   |     
        |                        |                  |will be ignored.                                          |
        --------------------------------------------------------------------------------------------------------
        |--cna_file              | REQUIRED         |Path to a tab-separated file containing copy number       |
        |                        |                  |aberrations (CNAs) where the first column of each line is |
        |                        |                  |a sample ID and subsequent columns contain gene names     |
        |                        |                  |followed by "(A)" or "(D)" indicating an amplification or |
        |                        |                  |deletion in that gene for the sample. Lines starting with |
        |                        |                  |'#' will be ignored.                                      |
        --------------------------------------------------------------------------------------------------------
        |--sample_file           | None             |Path to file listing sample IDs, one per line. Any SNVs or|
        |                        |                  |CNAs in samples not listed in this file will be ignored.  |
        |                        |                  |If HotNet is run with mutation permutation testing, all   |
        |                        |                  |samples in this file will be eligible for random mutations|
        |                        |                  |even if the sample did not have any mutations in the real |
        |                        |                  |data. If not provided, the set of samples is assumed to be|
        |                        |                  |all samples that are provided in the SNV or CNA data.     |
        --------------------------------------------------------------------------------------------------------
        |--gene_file             | None             |Path to file listing gene names, one per line. Mutations  |
        |                        |                  |in genes not listed in this file will be ignored. If      |
        |                        |                  |HotNet is run with mutation permutation testing, every    |
        |                        |                  |gene in this file will be eligible for random mutations   |
        |                        |                  |even if the gene did not have mutations in any samples in |
        |                        |                  |the original data. If not provided, the set of tested     |
        |                        |                  |genes is assumed to be all genes that have mutations in   |
        |                        |                  |either the SNV or CNA data.                               |
        --------------------------------------------------------------------------------------------------------
        |--min_freq              | 1                |The minimum number of samples in which a gene must have an|
        |                        |                  |SNV in order to be considered mutated in the heat score   |
        |                        |                  |calculation.                                              |
        --------------------------------------------------------------------------------------------------------
        |--cna_filter_threshold  | None             |Proportion of CNAs in a gene across samples that must     |
        |                        |                  |share the same CNA type in order for the CNAs to be       |
        |                        |                  |included. This must either be > .5, or the default, None, |
        |                        |                  |in which case all CNAs will be included.                  |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------


3. ###Delta selection###

    This step uses random data to select thresholds that should be used for edge weight removal in
    the HotNet run of step 4. The output of this step includes a list of 5 recommended deltas;
    users are advised to run HotNet with each of the 5 recommended deltas and manually select the
    best results.

    The random data can be either:

    * Permuted heat scores (for arbitrary data types),  

      or

    * Permuted mutation data (recommended for the case when mutation data is used to generate heat
      scores)
    
    The Python script `findThreshold.py` can be used to run the delta selection procedure. The
    required and optional parameters to the script are described below.

    For the permuted heat scores test, heat scores will simply be shuffled among genes. In this
    case, the first parameter to `findTheshold.py` should be `heat`, e.g.:

            python findThreshold.py heat <additional_parameters>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file       | REQUIRED         |Path to .mat file containing influence matrix.            |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | Li               |Variable name of the influence matrx in the .mat file.    |
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Path to tab-separated file containing an index in the     |
        |                        |                  |first column and the name of the gene represented at that |
        |                        |                  |index in the second column of each line.                  |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |JSON heat score file generated via generateHeat.py        |
        --------------------------------------------------------------------------------------------------------
        |-ms/--min_heat_score    | See description  |Minimum heat score for a gene to be eligible for inclusion|
        |                        |                  |in a returned connected component in permuted datasets. By|
        |                        |                  |default, all genes with positive heat scores will be      |
        |                        |                  |included. Note that if heat permutation is used to        |
        |                        |                  |generate permuted data sets, heat scores will be permuted |
        |                        |                  |among all genes that have scores, even those with scores  |
        |                        |                  |below the minimum for inclusion in output CCs.            |
        --------------------------------------------------------------------------------------------------------
        |-pgf/                   | None             |Path to file containing a list of additional genes that   |
        |--permutation_genes_file|                  |can have permuted heat values assigned to them in         |
        |                        |                  |permutation tests.                                        |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------

    For the permuted mutation data test, SNVs will be randomly generated in genes according to a
    specified background mutation rate, and each CNA block will be randomly placed on an equally-
    sized set of genes in the same chromosome. In this case, the first parameter to `findTheshold.py`
    should be `mutations`, e.g.:

        	python findThreshold.py mutations <additional_parameters>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        -------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file       | REQUIRED         |Path to .mat file containing influence matrix.            |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | Li               |Variable name of the influence matrx in the .mat file.    |
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Path to tab-separated file containing an index in the     |
        |                        |                  |first column and the name of the gene represented at that |
        |                        |                  |index in the second column of each line.                  |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |JSON heat score file generated via generateHeat.py        |
        --------------------------------------------------------------------------------------------------------
        |-glf/--gene_length_file | REQUIRED         |Path to tab-separated file containing gene names in the   |
        |                        |                  |first column and the length of the gene in base pairs in  |
        |                        |                  |the second column                                         |
        --------------------------------------------------------------------------------------------------------
        |-gof/--gene_order_file  | REQUIRED         |Path to file containing tab-separated lists of genes on   |
        |                        |                  |each chromosme, in order of their position on the         |
        |                        |                  |chromosome, one chromosome per line                       |
        --------------------------------------------------------------------------------------------------------
        |-b/--bmr                | REQUIRED         |Default background mutation rate                          |
        --------------------------------------------------------------------------------------------------------
        |-bf/--bmr_file          | None             |File listing gene-specific BMRs. If none, the default BMR |
        |                        |                  |will be used for all genes.                               |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------

4. ###HotNet run###

    This step performs the core HotNet algorithm: calculating a weighted graph based on influence
    matrix and heat score, removing edges with weight less than delta, and extracting the resulting
    connected components.

    The Python script `runHotNet.py` can be used to perform the HotNet run.  The required and
    optional parameters to the script are described in the tables below.  In addition to performing
    the run of the core algorithm, `runHotnet.py` also runs the significance testing described in
    step 5.

    To run HotNet without statistical significance testing, the final parameter to `runHotNet.py`
    should be `none`, e.g.:

            python runHotNet.py <additional_parameters> none

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file       | REQUIRED         |Path to .mat file containing influence matrix.            |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | Li               |Variable name of the influence matrices in the .mat files.|
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Path to tab-separated file containing an index in the     |
        |                        |                  |first column and the name of the gene represented at that |
        |                        |                  |index in the second column of each line.                  |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |JSON heat score file generated via generateHeat.py        |
        --------------------------------------------------------------------------------------------------------
        |-ms/--min_heat_score    | See description  |Minimum heat score for a gene to be eligible for inclusion|
        |                        |                  |in a returned connected component. By default, all genes  |
        |                        |                  |with positive heat scores will be included. Note that if  |
        |                        |                  |heat permutation is used to generate permuted data sets   |
        |                        |                  |for significance testing, heat scores will be permuted    |
        |                        |                  |among all genes that have scores, even those with scores  |
        |                        |                  |below the minimum for inclusion in output CCs.            |
        --------------------------------------------------------------------------------------------------------
        |-d/--delta              | REQUIRED         |Weight threshold for edge removal.                        |
        --------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size      | 3                |Minimum size connected components that should be returned.|
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | REQUIRED         |Output directory. Files results.json, components.txt, and |
        |                        |                  |significance.txt will be generated                        |
        --------------------------------------------------------------------------------------------------------

5. ###Statistical significance testing###

    This step calculates the statistical significance of obtaining the observed number of connected
    components of various sizes by comparing the observed number of connected components to the
    expected number from permuted data. As with delta selection, either heat scores or mutation data
    can be be permuted to generate the random data sets.

    As mentioned above, statistical significance testing is included as part of `runHotNet.py`. The
    additional required and optional parameters for including significance testing in the HotNet
    run are described below.

    To run statistical significance testing with permuted heat scores, include `heat` as a parameter
    after the parameters described in step 4, then include any additional parameters for the
    permutation test, e.g.:

            python runHotNet.py <params_from_step_4> heat <additional_significance_testing_params>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-pgf/                   | None             |Path to file containing a list of additional genes that   |
        |--permutation_genes_file|                  |can have permuted heat values assigned to them in         |
        |                        |                  |permutation tests.                                        |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------

    To run statistical significance testing with permuted mutation data, include `mutations` as a
    parameter after the parameters described in step 4, then include any additional parameters for
    the permutation test, e.g.:

            python runHotNet.py <params_from_step_4> mutations <additional_significance_testing_params>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-glf/--gene_length_file | REQUIRED         |Path to tab-separated file containing gene names in the   |
        |                        |                  |first column and the length of the gene in base pairs in  |
        |                        |                  |the second column                                         |
        --------------------------------------------------------------------------------------------------------
        |-gof/--gene_order_file  | REQUIRED         |Path to file containing tab-separated lists of genes on   |
        |                        |                  |each chromosme, in order of their position on the         |
        |                        |                  |chromosome, one chromosome per line                       |
        --------------------------------------------------------------------------------------------------------
        |-b/--bmr                | REQUIRED         |Default background mutation rate                          |
        --------------------------------------------------------------------------------------------------------
        |-bf/--bmr_file          | None             |File listing gene-specific BMRs. If none, the default BMR |
        |                        |                  |will be used for all genes.                               |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------

6. ###Visualization###

    You can visualize the subnetworks output by HotNet using the `makeResultsWebsite.py` script.
    It takes the following parameters:

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--results_files      | REQUIRED         |Paths to results.json files output by HotNet. Multiple    |
        |                        |                  |file paths may be passed.                                 |
        --------------------------------------------------------------------------------------------------------
        |-ef/--edge_file         | REQUIRED         |Path to TSV file listing edges of the interaction network,|
        |                        |                  |where each row contains the indices of two genes that are |
        |                        |                  |connected in the network.                                 |
        --------------------------------------------------------------------------------------------------------
        |-nn/--network_name      | Network          |Display name for the interaction network.                 |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | REQUIRED         |Output directory.                                         |
        -------------------------------------------------------------------------------------------------------- 

    To view the resulting visualizations, navigate to the output directory and run
    `python -m SimpleHTTPServer`, then visit `http://localhost:8000` in a browser.

7. ###Output annotation###

    The Python script `annotateOutput.py` marks each gene in the output subnetworks with the number
    of samples that contained a mutation in the gene.  It also includes the total number and
    percentage of samples mutated in each subnetwork.  Note that this output annotation is only
    applicable for HotNet runs using mutation data.  The following parameter is required:

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |--hotnet_output_json    | REQUIRED         |Path to JSON output file produced by running runHotNet.py |
        --------------------------------------------------------------------------------------------------------

    The output from the script will be generated as `annotated_components.txt` in the same output
    directory specified for `runHotNet.py`.


Passing parameters
------------------------
To avoid the need to continually pass a large number of parameters on the command line, some or all
parameters to the scripts described above can be specified via a configuration file by passing
`@<ConfigFileName>` as a command-line parameter, e.g. `python runHotNet.py @testConf.txt`.  The
configuration file simply contains additional parameters as they would be specified on the
command line (though, unlike on the command line, line breaks are permitted between parameters in
config files). Multiple configuration files can be used, and configuration files can be combined
with parameters specified directly on the command line. If the same parameter is specified in
multiple configuration files or in both a config file and on the command line, the last set value
will be used.


Influence matrices
------------------------
The following pre-computed influence matrices are provided:

* matrices for HPRD v9 (http://www.hprd.org/)
* matrices for iRefIndex v9 (http://irefindex.uio.no/wiki/iRefIndex)

Due to their large size, the influence matrices are not distributed with the code but are available at:
http://compbio.cs.brown.edu/projects/hotnet/

Sample data
------------------------
In the 'example' folder we provide sample data for:

* snv_file: example.snv
* cna_file: example.cna
* heat_file:  example.heat
* gene\_file: tested\_genes.txt
* sample\_file: tested\_samples.txt

We also include in the 'auxiliary_files' folder some of the files that can be used with HotNet
when mutation data from whole exome/genome sequencing and from genome-wide copy number assays is
analyzed.  NOTE: these files are not intended for use when targeted sequencing data is analyzed,
since they assume that all genes have been tested for mutations and copy number aberrations.

Files:

* order\_genes\_hg19\_gistic\_bychrom.txt: order of genes in chromosomes, for `-gof/--gene_order_file`
    parameter in runHotNet.py
* all\_human\_genes.txt: list of human genes, to be used as list of tested genes for `--gene_file`
    parameter in generateHeat.py. Note that genes (e.g., TTN) that cause known sequencing/mapping
    artifacts have been removed from the list.
* gene\_length\_2013-08-09.txt: lengths of genes (taken from Ensembl on 2013-08-09) for
    `-glf/--gene_length_file` in findThreshold.py and runHotNet.py.

Finally, the 'configs' directory contains configuration files that can be used to run HotNet on the
sample data using HPRD.  Note that these configuration files assume that the HPRD influence matrix
and gene-index files have been downloaded and saved with their default names in the 'influence_matrices'
directory.  The samples can be run as follows:

        python generateHeat.py @example/configs/heat.config
        python findThreshold.py @example/configs/delta.config
        python runHotNet.py @example/configs/run.config
        python runHotNet.py @example/configs/significance.config

Citation
------------------------
If you use HotNet in your work, please cite:

F. Vandin, E. Upfal, and B.J. Raphael. (2011) Algorithms for Detecting Significantly Mutated
Pathways in Cancer. Journal of Computational Biology. 18(3):507-22.

F. Vandin, P. Clay, E. Upfal, and B. J. Raphael. Discovery of Mutated Subnetworks Associated with Clinical Data in Cancer. In Proc. Pacific Symposium on Biocomputing (PSB), 2012.
