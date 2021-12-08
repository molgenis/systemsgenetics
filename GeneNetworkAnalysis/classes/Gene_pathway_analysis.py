"""
File:         Gene_pathway_analysis.py
Created:      2020-04-29
Last Changed: 2021-01-08
Author(s):    H.H.Wiersma

Code to create distribution plot and a mask image used in the research plan
Copyright (C) 2019 H.H.Wiersma

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

from utilities import timer_print, read_pd_df, logging_print, \
    create_auc_output_file, remove_gene_from_pathway

import os
import glob
import time
from datetime import datetime
import logging
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import scipy.special
import multiprocessing as mp
import queue
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression

analysis_types = {
    "REGRESSION": "regression",
    "T_TEST": "t_test",
}

class Gene_Pathway_Analysis:
    # Class to perform the pathway analysis which calculate the gene
    # contribution per pathway based on a t-test approach (old gene network
    # implementation) or based on a logistic model fit.
    def __init__(self, components_data_path, matrix_path,
                 output_dir, analysis_type, minimal_number_of_genes=None,
                 background_genes_path=None, n_cores=None,
                 split_start=None, split_end=None, multi_node_node_id=None,
                 multi_node_num_nodes=None, multi_node_output_dir=None,
                 force=False, output_disable_txt=False,
                 output_disable_gzip=False, output_disable_pickle=False
                 ):
        """
        Class to perform the pathway analysis
        :param components_data_path: The path of the input file. The input file contains the gene expression values where the rows correspond to the samples and the columns to the genes. The input file can be given as pandas pickled dataframe or as tab separated txt file. The txt input file is automatically decompressed if the file paths ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’.
        :param matrix_path: The path to the pathway matrix file(s). The parameter must contain at least one path to a pathway matrix but can also contain multiple matrixes by separating the paths by commas. <br> The input file(s) can be given as pandas pickled dataframe(s) or as tab separated txt file(s). The txt input file is automatically decompressed if the file path ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’
        :param output_dir: Path to the output directory. The output files will be created in this directory. The directory will automatically be created if not exists
        :param analysis_type: Set the analysis which will be used to fit the model and calculate the gene pathway scores. Available analysis types: “REGRESSION” and “T_TEST”. <br><br> REGRESSION: The updated method which used the logistic regression model to fit an model which is used to calculate the gene z_scores by using permutation of the input matrix. <br><br> T_TEST: The original gene network method which use an t-test apparoach combined with an linear model to calculate the gene p-values which are translated to the gene z-scores
        :param minimal_number_of_genes: The minimal number of genes which must be present in the pathway to be analyzed. Pathways with a lower number of genes are skipped
        :param background_genes_path: The path to the background gene file(s). The genes mentioned in the file(s) are used to fit the model instead of all the genes present in the input matrix. If multiple pathway matrixes were given, than the background files must be given in the same order as the pathway matrix. <br> The parameter must contain at least one path to a pathway matrix but can also contain multiple matrixes by separating the paths by commas. <br> The input file(s) can be given as pandas pickled dataframe(s) or as tab separated txt file(s). The txt input file is automatically decompressed if the file path ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’
        :param n_cores: Number of cores which will be used to process the pathways on different processor cores. Default is the number of available
        :param split_start: The first index of the pathway matrix to start the analysis. Can be used to subset the input pathway matrix (together with --end)
        :param split_end: The last index of the pathway matrix to end the analysis. Can be used to subset the input pathway matrix (together with --start)
        :param multi_node_node_id: The id of the single node (zero first). This option can be used to execute the script in an multi node environment (cluster). Every node will execute an equal number of pathways. Option must be used together with the parameter “num_nodes”. <br><br> For example, if the cluster has 10 nodes and the pathway file contains 100 pathways, every node will process 10 pathways. The first node (--node 0) will process the first 10 pathways, the second node (--node 1) will process pathway 11 until 20, the third node (--node 2) will process pathway 21 until 30, etc
        :param multi_node_num_nodes: Total number of nodes used in the cluster (multi node) environment. This option can be used to execute the script in an multi node environment (cluster). Every node will execute an equal number of pathways.
        :param multi_node_output_dir: The overall output dir of the analysis. If the analysis is excecuted in a multi node environment, this path will be the head direcory of the analysis. The output directory variable are the node specific output directory
        :param force: Do not use the cached versions of the input data but the original input files
        :param output_disable_txt: Do not export the output as text files
        :param output_disable_gzip: Do not gzip the data output text files
        :param output_disable_pickle: Do not export the output as pickle Pandas files
        """

        # Set the parameters
        self.start_time = time.time()
        self.components_data_path = components_data_path
        self.output_dir = output_dir
        self.analysis_type = analysis_type
        self.matrix_path = matrix_path
        self.background_genes_path = background_genes_path

        # The minimal number of genes in a gene set must be three, so set
        # three if the no number of genes is given
        self.minimal_number_of_genes = -1
        if minimal_number_of_genes is not None and minimal_number_of_genes != "":
            self.minimal_number_of_genes = int(minimal_number_of_genes)
        if self.minimal_number_of_genes < 3:
            self.minimal_number_of_genes = 3
            logging_print("Minimal number of genes must be 3 of higher, so the value is set to 3")

        # The paramters below are used for the multiprocessing
        self.split_start = split_start
        self.split_end = split_end
        self.multi_node_node_id = multi_node_node_id
        self.multi_node_num_nodes = multi_node_num_nodes
        self.multi_node_output_dir = multi_node_output_dir

        self.n_cores = mp.cpu_count()
        if n_cores is not None and n_cores != "":
            self.n_cores = int(n_cores)

        self.components_data = None
        self.background_genes_data = None
        self.matrix_data = None
        self.components_data_permutated = None

        # results
        self.gene_pathway_count = None
        self.pathway_gene_scores = None
        self.auc_values = None
        self.p_values = None
        self.bonf_p_values = None
        self.reject = None
        self.force=force
        self.method_stats = {}

        # output
        self.output_disable_txt = output_disable_txt
        self.output_disable_gzip = output_disable_gzip
        self.output_disable_pickle = output_disable_pickle

    def run(self):
        ###
        ### Method to excecute the complete gene-pathway analysis
        ###

        # Write the basic information about the analysis to the log file
        self.write_base_loginfo()
        # Read the input files
        self.read_files()
        # Perform gene intersection between all the input files to match the
        # input dataframes for further analyse
        self.perform_gene_intersection()
        # Remove all the pathways with less than the given number of
        # gene annotations
        self.perform_filtering()
        # Split up the input data for multiprocessing
        self.perform_multinode_processing()
        # Caluclate the number of gene annotations per pathway
        self.calculate_gene_pathway_count()
        # Create or load the gene/component permutation matrix
        self.create_permutations_matrixes()
        # Perform the pathway analysis
        self.perform_analysis()
        # Calculate the pathway AUC values
        self.calculate_auc_values()
        # Calculate the pathway p-value
        self.calculate_wilcox_p_value()

        # Calculate the bonferroni corrected p-values of the analyse data if
        # the analyse was not executed with multi node processing.
        # If multi processing is used, the bonferroni corrected p-values will
        # be calculated over te overall results of all the nodes.
        if self.multi_node_num_nodes == None :
            self.calculate_p_value_bonferroni_correction()

        # Create the output files
        self.write_output_files()
        # Remove all the temp files
        self.remove_temp_files()

        # Merge the output files from the different nodes and calculate the
        # bonferroni corrected pathway p-values over all the analysed pathways
        if self.multi_node_num_nodes != None:
            self.merge_output_files()

        # Write the final information to the log file
        self.write_last_loginfo()

    def write_base_loginfo(self):
        ###
        ### Write basic log info of the analysis to the logfile
        ###
        info = "## START GENE PATHWAY ANALYSIS (new p value method) ##\n" \
               "DATE:\t{date}\n" \
               "Component file:\t{component_file}\n" \
               "Matrix file:\t{matrix_file}\n" \
               "Background gene file:\t{background_gene_file}\n" \
               "Output dir:\t{output_dir}\n" \
               "Analysis type:\t{analysis_type}\n" \
               "Number of cores:\t{num_cores}\n" \
               "".format(date=datetime.now(),
                         component_file=self.components_data_path,
                         matrix_file=self.matrix_path,
                         background_gene_file=self.background_genes_path,
                         output_dir=self.output_dir,
                         analysis_type=self.analysis_type,
                         num_cores=self.n_cores)
        logging_print(info)

    def write_last_loginfo(self):
        ###
        ### Write last line of the logfile
        ###
        timer_print(self.start_time, prefix="## ANALYSIS READY")

    def read_files(self):
        ###
        ### Method to read the input files
        ###
        rf_start_time = time.time()

        # read input component file
        if os.path.isfile(self.components_data_path):
            self.components_data = read_pd_df(self.components_data_path,
                                              {
                                                  "sep": "\t",
                                                  "index_col": 0
                                              },
                                              force=self.force)
            components_data_info = "Components dataframe n_row: {}, " \
                             "n_col: {}\n" \
                             "first column headers: {}\n" \
                             "first row index: {}".format(
                *self.components_data.shape,
                ", ".join(
                    self.components_data.columns.values[
                    :5]),
                ", ".join(
                    self.components_data.index.values[
                    :5]),
                )
            logging_print(components_data_info)
        else:
            raise FileNotFoundError("Cannot find input file: {}".format(
                self.components_data_path
            ))

        # read matrix file
        self.matrix_data = read_pd_df(self.matrix_path,
                                      {
                                          "sep": "\t",
                                          "index_col": 0
                                      },
                                      force=self.force,
                                      proc_df_before_save=lambda df: df == 1.0
                                      )

        matrix_data_info = "Matrix dataframe n_row: {}, " \
                         "n_col: {}\n" \
                         "first column headers: {}\n" \
                         "first row index: {}".format(
            *self.matrix_data.shape,
            ", ".join(
                self.matrix_data.columns.values[
                :5]),
            ", ".join(
                self.matrix_data.index.values[
                :5]),
            )
        logging_print(matrix_data_info)

        # reading background gene path
        if self.background_genes_path is not None and self.background_genes_path != '':
            self.background_genes_data = read_pd_df(self.background_genes_path,
                                                    {
                                                        "sep": "\t",
                                                        "header": None
                                                    },
                                                    force=self.force
                                                    )

            # convert to serie instead of matrix
            self.background_genes_data = self.background_genes_data.iloc[:, 0]

            background_genes_data_info = "Background genes file loaded:\n" \
                               "number of genes: {}, " \
                               "first genes: {}\n".format(
                self.background_genes_data.shape[0],
                ", ".join(
                    self.background_genes_data.values[
                    :5])
            )
            logging_print(background_genes_data_info)

        timer_print(rf_start_time, prefix="Reading input file ready")

    def has_background_genes(self):
        ##
        ## Method to check if background file is loaded
        ##

        return self.background_genes_data is not None

    def perform_gene_intersection(self):
        ###
        ### Perform gene intersection so that alle
        # the matrixes contains the same genes
        ###

        pgi_start_time = time.time()
        intersect_genes = self.components_data.index.intersection(self.matrix_data.index)
        self.components_data = self.components_data.loc[intersect_genes, :]
        self.matrix_data = self.matrix_data.loc[intersect_genes, :]
        if self.has_background_genes():
            self.background_genes_data = pd.Series(list(set(self.background_genes_data).intersection(set(intersect_genes))))

        timer_print(pgi_start_time, prefix="Gene intersection ready")

    def perform_filtering(self):
        ###
        ### Remove all the pathways which contains less genes
        # than the minimal number of genes set by the parameters
        ###

        if self.minimal_number_of_genes > 0:
            pf_start_time = time.time()
            matrix_selection = self.matrix_data.sum(axis=0) >= self.minimal_number_of_genes
            self.matrix_data  = self.matrix_data.loc[:, matrix_selection[matrix_selection].index]
            logging_print("Minimal gene filtering: {} number of pathways over".format(self.matrix_data.shape[1]))
            timer_print(pf_start_time, prefix="Minimal gene in pathway filtering ready")

    def create_permutations_matrixes(self):
        ###
        ### Mehtod to create the permutation matrix of the orginal
        ### components data frame which will be used to calculate the z-scores
        ### this matrix will be used after the model fit to center and scale
        ### the gene logodds score to a z score
        ###
        ### To guarantee the same results during different runs,
        # the permutation matrix will be saved next to the components matrix
        # file with the file name {component_file_name}_permutation_matrix.pkl.
        # If this file is already present, than the file is loaded instead
        # of creating a new permutation matrix. This loading process can be
        # disabled by setting the force parameter to true, than a
        # new matrix is created
        #
        # The permutation metrix will be created by random suffeling
        # of the index of the compoenents dataframe
        ###

        perm_start_time = time.time()
        permutation_path = "{}_permutation_matrix.pkl".format(self.components_data_path)

        # check if permutation matrix is not already present
        if os.path.isfile(permutation_path) and self.force == False:
            self.components_data_permutated = pd.read_pickle(permutation_path)
            logging_print("Permutation dataframe is loaded from '{}', "
                          "shape: {}".format(permutation_path,
                                             self.components_data_permutated.shape))
        else:
            # create the permutation matrix by random sort of the
            # index of the dataframe. This is done by the permutation method of numpy
            self.components_data_permutated = self.components_data.set_index(
                np.random.permutation(self.components_data.index))
            self.components_data_permutated.to_pickle(permutation_path)
            logging_print("Permutation table is created, shape: {}".format(
                self.components_data_permutated.shape))
        timer_print(perm_start_time, prefix="permutation is ready")

    def perform_multinode_processing(self):
        ###
        ### The script contains functionality to split up the data over
        # multiple nodes to speedup the analysis process.
        # This means that every node process only a small number of pathways
        # from the pathway matrix. This method selects the pathways which
        # will be processed by the current node
        ###

        # First subset the data by the user parameters
        if self.split_end is not None:
            self.matrix_data = self.matrix_data.iloc[:, :self.split_end]

        if self.split_start is not None:
            self.matrix_data = self.matrix_data.iloc[:, self.split_start:]

        # log this subset information
        if self.split_start is not None or self.split_end is not None:
            logging_print("Trim matrix, new start: {start}, "
                       "new end: {end}, dataframe size: {df_size}".format(
                start=self.split_start,
                end=self.split_end,
                df_size=self.matrix_data.shape
            ))

        # handle multi node processing on a cluster
        if self.multi_node_num_nodes is not None and \
                self.multi_node_node_id is not None:
            # Calulate the number of pathways which must be processed per node
            num_pathways = self.matrix_data.shape[1]
            pathways_per_node = num_pathways // self.multi_node_num_nodes
            if pathways_per_node * self.multi_node_num_nodes < num_pathways:
                pathways_per_node += 1

            # calculate the start and stop positions of the pathway matrix
            # based on the number of pathways processed
            # per node and the node id given by as parameter
            start_id = self.multi_node_node_id * pathways_per_node
            end_id = (self.multi_node_node_id + 1) * pathways_per_node
            if end_id > num_pathways:
                end_id = num_pathways

            # Subset the matrix for the multinode processing
            self.matrix_data = self.matrix_data.iloc[:, start_id:end_id]

            # log this multi node subset information
            logging_print("## Process on multiple nodes ##\n"
                          "Node: {node_id} of {node_num}\n"
                          "Number pathways per node: {num_pathways_per_node} of {num_pathways}\n"
                          "Start id: {start_id}\n"
                          "End id: {end_id}\n"
                          "Dataframe size: {df_size}".format(
                node_id=self.multi_node_node_id + 1,
                node_num=self.multi_node_num_nodes,
                num_pathways_per_node=pathways_per_node,
                num_pathways=num_pathways,
                start_id=start_id,
                end_id=end_id,
                df_size=self.matrix_data.shape
            ))

    def calculate_gene_pathway_count(self):
        ###
        ### Return the number of genes per pathway
        ###

        self.gene_pathway_count = self.matrix_data.sum(axis=0)

    def perform_analysis(self):
        ###
        ### This Method performs starts the analysis per pathway.
        # The method submit the single pathway jobs to multiple workers
        # because the work is splitted up over the different
        # cores of the machine. The method used a queue with which contains
        # a "job" for every pathway. The Takes the a job from the queue,
        # processed it and push the results in a results queue, and starts
        # again from the start by getting a new job from the queue.
        # This is repeated until all the jobs are processed from the queue.
        # The worker puts a "DONE" string in the result queue if the
        # worker is done and there are no more jobs in the input queue.
        #
        # the new results will be saved every two minutes in a temporary file
        # to prevent result loss after crash. The temp file will be saved
        # in the output directory with the filename:
        # temp_results_analysis_z_scores_{date}.pkl and will
        # be automatically loaded when the program is (re)started again.
        #
        # after the analysis, the results will also be temporary stored
        # in the temp_pathway_gene_scores_temp_{date}.pkl file
        #
        ###

        pa_start_time = time.time()
        # results from the analysis
        total_z_score_results = []

        # check if some temp files where present
        temp_z_score_paths = glob.glob(os.path.join(self.output_dir,
                                                    "temp_results_analysis_z_scores_*.pkl"))
        temp_z_scores = None
        temp_already_processed_pathways = None
        if len(temp_z_score_paths) > 0:
            # load the data from the temporary file
            temp_z_scores = pd.read_pickle(temp_z_score_paths[0])
            logging_print("Temp file '{}' with already processed pathways loaded. size df: {}".format(temp_z_score_paths[0], temp_z_scores.shape))
            temp_already_processed_pathways = list(temp_z_scores.columns.to_numpy())

        # manager and queue for the input data
        pathway_manager = mp.Manager()
        pathway_queue = pathway_manager.Queue()

        # manager and queue for the results
        retults_manager = mp.Manager()
        results_queue = retults_manager.Queue()

        # Add the Pathways tot the input queue if the pathway
        # was not already processed and saved in the temp file
        for index, row in self.matrix_data.iteritems():
            if temp_z_scores is not None:
                # Check the temp file
                if index in temp_already_processed_pathways:
                    # Put the results from the temp file in the results list
                    total_z_score_results.append(temp_z_scores.loc[:, index])
                else:
                    # Pathway is not already processed, so put the
                    # pathway index in the queue
                    pathway_queue.put(index)
            else:
                pathway_queue.put(index)
        # log the number of already processed pathways
        logging_print("total pathways already done: {}".format(
            len(total_z_score_results)))

        # Setup the workers for every core. One core will not be used
        # for processing but will be used for head job
        n_workers = self.n_cores -1
        processes = []
        for _ in range(n_workers):
            processes.append(mp.Process(target=single_pathway_worker,
                                        args=(self.components_data,
                                              self.matrix_data,
                                              pathway_queue,
                                              results_queue,
                                              self.background_genes_data,
                                              self.analysis_type,
                                              -1,
                                              self.components_data_permutated
                                              )))
        # start the workers
        for process in processes:
            process.start()

        #
        # Process the results from the workers
        #
        total_done = 0
        last_save_time = time.time()
        save_time_in_minutes = 2
        z_score_file_path = None
        # Do this until all the workers are ready
        while True:
            try:
                if total_done >= n_workers:
                    # All workers are ready
                    break
                if results_queue.empty():
                    # Wait for the next results
                    time.sleep(5)
                else:
                    # Process the results

                    # get the result from the result queue
                    results = results_queue.get(True, timeout=1)
                    if isinstance(results, str) and results == "DONE":
                        # One worker is done
                        total_done += 1
                        logging_print("total workers done {} of {}".format(total_done, n_workers))
                    else:
                        # Save the results
                        total_z_score_results.append(results)

                    # save temp results
                    if last_save_time + 60 * save_time_in_minutes < time.time():
                        logging_print("save temp results: {}".format(datetime.now()))
                        last_save_time = time.time()
                        temp_dataframe_z_scores = pd.DataFrame(total_z_score_results).T
                        new_z_score_file_path = os.path.join(self.output_dir, "temp_results_analysis_z_scores_{}.pkl".format(datetime.now()))

                        # create the new log file
                        temp_dataframe_z_scores.to_pickle(new_z_score_file_path)

                        # remove an earlier logfile
                        if z_score_file_path is not None and os.path.isfile(z_score_file_path):
                            os.remove(z_score_file_path)
                        z_score_file_path = new_z_score_file_path

            except queue.Empty:
                time.sleep(1)
                continue
        # all workers are done, remove the workers
        for process in processes:
            process.join()

        # save the results in a new temp file
        print("all processes ready")
        self.pathway_gene_scores = pd.DataFrame(total_z_score_results).T

        pathway_gene_scores_temp_file_path = os.path.join(self.output_dir,
                                             "temp_pathway_gene_scores_temp_{}.pkl".format(
                                                 datetime.now()))
        self.pathway_gene_scores.to_pickle(pathway_gene_scores_temp_file_path)

        # log the total processing time
        timer_print(pa_start_time, prefix="Performing analysis ready")

    def calculate_auc_values(self):
        ###
        ### Method to claculate the AUC value of the pathway predictions.
        # The AUC is calculated by using the ROC_curve value between
        # the calculated pathway z scores and the original gene pathway values.
        # From the COV curve the AUC value is calculated
        ###

        cav_st_stime = time.time()
        auc_values = {}
        # Calculate the AUC for every pathway
        for pathway_id in self.pathway_gene_scores.columns:
            if pathway_id in list(self.matrix_data.columns):
                # get the datasets
                pathway_gene_set = self.matrix_data.loc[:, pathway_id]
                pathway_gene_scores = self.pathway_gene_scores.loc[self.matrix_data.index, pathway_id]
                try:
                    # calculate the ROC curve
                    log_reg_fpr, log_reg_tpr, _ = roc_curve(pathway_gene_set * 1,
                                                            pathway_gene_scores)
                    # calculate the corresponded AUC value
                    auc_values[pathway_id] = auc(log_reg_fpr, log_reg_tpr)
                except ValueError:
                    logging_print("Value error in auc calculation of pathway: {}".format(pathway_id))
            else:
                logging_print("HPO term: {} not in matrix".format(pathway_id))

        # safe the AUC values
        auc_calc_values = pd.Series(auc_values)

        # replace AUC values higher dan 1 by 1
        auc_calc_values[auc_calc_values > 1.0] = 1.0
        self.auc_values = auc_calc_values
        timer_print(cav_st_stime, prefix="AUC value calculation ready")

    def calculate_wilcox_p_value(self):
        ###
        ### Calulate the P value corresponend to the Pathway prediction and
        # the AUC value of the prediction. The Wilcox test is used to
        # calculate the p value between the predicted gene values of the
        # genes which is present in the pathway given by the original
        # gene pathway matrix against the predicted gene pathway scores
        # which are not present in the pathway given by the oroginal
        # gene pathway matrix
        ###

        cwpv_st_time = time.time()
        wal_p_values = {}

        # calculate the p-value for every pathway
        for pathway_id in self.auc_values.index:
            if pathway_id in list(self.matrix_data.columns):
                # get the gene pathway information from the original given matrix
                selected_pathway = self.pathway_gene_scores.loc[self.matrix_data.index, pathway_id]

                # get the gene prediction scores of the genes which are
                # included in the pathway given by the original pathway matrix
                include_in_pathway = selected_pathway[self.matrix_data.loc[:, pathway_id]]

                # get the gene prediction scores of the genes which are not
                # included in the pathway given by the original pathway matrix
                include_not_in_pathway = selected_pathway[~self.matrix_data.loc[:, pathway_id]]
                try:
                    # calculate the p value
                    _, p_value = scipy.stats.mannwhitneyu(include_in_pathway,
                                                       include_not_in_pathway,
                                                       use_continuity=True,
                                                       alternative="two-sided")
                    wal_p_values[pathway_id] = p_value
                except ValueError:
                    logging_print("Value error in p value calculation, pvalue 0.0 is set to pathway: {}".format(pathway_id))
                    wal_p_values[pathway_id] = 0.0

        # save the p-values
        self.p_values = pd.Series(wal_p_values)
        timer_print(cwpv_st_time, prefix="Wilcox p value calculation ready")

    def calculate_p_value_bonferroni_correction(self, alpha=0.05):
        ###
        ### Method to perform the Bonferonni correction of the
        # pathway prediction p values
        ###

        cpvbc_st_time = time.time()
        self.bonf_p_values = self.p_values * self.p_values.shape[0]
        timer_print(cpvbc_st_time,
                         prefix="Bonferroni p value correction ready")

    def write_output_files(self):
        ###
        ### Method to write the output files
        ###

        wof_start_time = time.time()

        # Create AUC output file
        create_auc_output_file(self.auc_values,
                               self.gene_pathway_count,
                               self.p_values,
                               self.output_dir,
                               self.bonf_p_values,
                               txt_gzip=self.output_disable_gzip == False,
                               export_txt=self.output_disable_txt == False,
                               export_pkl=self.output_disable_pickle == False
                                          or self.multi_node_num_nodes != None
                               )

        # Export gene pathway scores
        if self.output_disable_txt == False:
            if self.output_disable_gzip:
                self.pathway_gene_scores.to_csv(
                    os.path.join(self.output_dir, "gene_pathway_scores.txt"),
                    sep='\t')
            else:
                self.pathway_gene_scores.to_csv(
                    os.path.join(self.output_dir, "gene_pathway_scores.txt.gz"),
                    sep='\t')
        if self.output_disable_pickle == False or \
                self.multi_node_num_nodes != None:
            self.pathway_gene_scores.to_pickle(
                os.path.join(self.output_dir, "gene_pathway_scores.pkl"))

        timer_print(wof_start_time, prefix="Writing output files ready")

    def remove_temp_files(self):
        ###
        ### Method to remove the exsiting temp files after analysis
        ###

        for file_path in glob.glob(os.path.join(self.output_dir,
                                                "temp_results_analysis_z_scores_*.pkl")):
            os.remove(file_path)

    def merge_output_files(self):
        ###
        ### Method to merge the output files from the different nodes.
        # If all the analysis excecuted on the different nodes are ready,
        # than all the input files will be read by this method and
        # merged together. The method also adds a bonferroni significant
        # p-value of the gene p-values to the output file which will be
        # corrected based on the total number of analysed pathways.
        ###

        mof_start_time = time.time()

        # Check if the script is excecuted over multiple nodes
        if self.multi_node_output_dir:

            # check if all the output files of the nodes are present,
            # e.q. all the nodes are ready
            node_output_auc_pred_files = glob.glob(os.path.join(
                self.multi_node_output_dir, "*", "predictions_auc.pkl"))
            if len(node_output_auc_pred_files) == self.multi_node_num_nodes:
                logging_print("Merge output files")

                ###
                ### merge AUC prediction files
                ###
                comp_df_auc_list = []
                for node_output_auc_file_path in node_output_auc_pred_files:
                    # read the file
                    file_auc_df = pd.read_pickle(node_output_auc_file_path)
                    # add dataframe to list
                    comp_df_auc_list.append(file_auc_df)
                #merge the data frames
                overall_auc_df = pd.concat(comp_df_auc_list)

                # Calculate the bonferroni corrected p-value of the pathway
                # prediction p-value
                bonf_values = overall_auc_df["pValue"] * overall_auc_df.shape[0]
                bonf_values[bonf_values > 1] = 1
                bonf_values[bonf_values < 0] = 0
                overall_auc_df["bonferroni"] = bonf_values

                # create the output text file
                if self.output_disable_txt == False:
                    output_file_path_csv = os.path.join(self.multi_node_output_dir,
                                                        "predictions_auc_bonf.txt.gz")
                    if self.output_disable_gzip:
                        output_file_path_csv = os.path.join(self.multi_node_output_dir,
                                                            "predictions_auc_bonf.txt")
                    overall_auc_df.to_csv(
                        output_file_path_csv,
                        columns=["geneCount", "pValue", "auc", "bonferroni"],
                        sep='\t')

                # create the output pickle file
                if self.output_disable_pickle == False:
                    overall_auc_df.to_pickle(os.path.join(
                        self.multi_node_output_dir, "predictions_auc_bonf.pkl"))

                ###
                ### merge gene pathway files
                ###
                node_output_gene_pathway_pred_files = glob.glob(
                    os.path.join(self.multi_node_output_dir, "*",
                                 "gene_pathway_scores.pkl"))
                comp_df_gene_pathway_list = []
                for node_output_auc_file_path in node_output_gene_pathway_pred_files:
                    # read the output file
                    file_gene_pathway_df = pd.read_pickle(node_output_auc_file_path)
                    comp_df_gene_pathway_list.append(file_gene_pathway_df)

                # Create the complete dataframe
                overall_gene_pathway_df = pd.concat(comp_df_gene_pathway_list,
                                                    axis=1)

                # create the output text file
                if self.output_disable_txt == False:
                    output_gene_pathway_pred_path = os.path.join(
                        self.multi_node_output_dir,
                        "gene_pathway_scores.txt.gz")

                    if self.output_disable_gzip:
                        output_gene_pathway_pred_path = os.path.join(
                            self.multi_node_output_dir,
                            "gene_pathway_scores.txt")
                    overall_gene_pathway_df.to_csv(
                        output_gene_pathway_pred_path,
                        sep='\t')

                # create the output pickle file
                if self.output_disable_pickle == False:
                    overall_gene_pathway_df.to_pickle(
                        os.path.join(self.multi_node_output_dir,
                                     "gene_pathway_scores.pkl"))
        timer_print(mof_start_time, prefix="Writing merged outputfile ready")


def single_pathway_worker(components_data, matrix_data, pathway_queue,
                          results_queue, background_genes_data=None,
                          analysis_type=None, max_end_time=-1,
                          components_data_permutated=None):
    ###
    ### This function is the worker which excecutes the basic functionality of
    # the worker. The worker will work until there are no more jobs to
    # process or as the maximal time is elapsed. The worker will get a
    # pathway id from the input queue, select the corresponded pathway
    # genes from the input table and will process that pathway to
    # calculate the gene pathway scores. These results will be pushed in the
    # results queue. When the results are saved, a new pathway id will be
    # extracted from the input queue and the process starts again. The input
    # queue will be empty if all the pathways are processed and the worker will
    # then return a string with "DONE" to let the head core known
    # that this worker is ready.
    ###

    try:
        while True:
            try:
                if pathway_queue.empty():
                    # no more pathways to process
                    break
                if 0 < max_end_time <= time.time():
                    # maximal end time is reached
                    break
                # get the pathway from the queue
                single_pathway_id = pathway_queue.get(True, timeout=1)

                # process the single pathway
                results = perform_single_pathway(components_data,
                                                 matrix_data.loc[:, single_pathway_id],
                                                 background_genes_data,
                                                 analysis_type,
                                                 components_data_permutated)
                # rename the result series and put the results
                # in the result queue
                results.name = single_pathway_id
                results_queue.put(results)
            except queue.Empty:
                time.sleep(1)
                continue

        # worker is done, put tht done string in the results queue
        results_queue.put("DONE")
        time.sleep(1)
        return

    except Exception as ex:
        # process the errors
        print("ERROR in worker: {}".format(ex))
        print("exception done")
        results_queue.put("DONE")
        time.sleep(1)
        logging.exception("ERROR in worker")


def single_row_analysis_t_test(row, components_data,
                               matrix_genes, background_genes_data,
                               overall_test, is_in_pathway=False):
    """
    Method to perform the t-test analysis of a single gene in a single pathway.
    :param row: Row of the gene pathway matrix
    :param components_data: Components data matrix
    :param matrix_genes: Boolean Serie which genes in active or not in the pathway
    :param background_genes_data: Serie with the genes names which must be used
    :param overall_test: The results of the overall t-test analysis between the component values of genes active against genes not active in the pathway
    :param is_in_pathway: Boolean if the gene is active in the original pathway description
    :return: correlation value, p value of the correlation

    If the gene is not active in the original description
    (is_in_pathway = false), the overall t-test
    results will be used instead of performing a new t-test. The Pearson
    correlation will be calculated between the original results and the
    components data and the results from the correlation are returned. If the
    gene is present in the original pathway description, the gene will be
    temporary set to false in the pathway description and a new t-test will
    be executed. The results from this t-test are translated to z-scores and
    the Pearson correlation is calculated between these z-scores and
    the components values.
    """

    try:
        # Check if the current gene to analyse is active in the original
        # pathway annotation
        if is_in_pathway:
            #gene is active

            # Remove the gene from the pathway (set the gene to false in the
            # annotation overview)
            new_gene_set = remove_gene_from_pathway(matrix_genes, row.name)
            # Perform gene intersection to match all the input Series
            intersect_components_data_keys = components_data.index.intersection(new_gene_set.index)

            # Perform a new t-test analysis
            z_scores, _, _ = perform_t_test(components_data.loc[intersect_components_data_keys, :],
                                            new_gene_set,
                                            background_genes=background_genes_data)
            if z_scores.sum() > 0:
                # Perform a correlation between the new calculated p-values
                # and the component values
                cor = scipy.stats.pearsonr(row, z_scores)
                return cor
        else:
            # gene is not in pathway
            if overall_test.sum() > 0:
                # Perform a correlation between the p-values results from the
                # overall t-test analysis and the component values
                cor = scipy.stats.pearsonr(row, overall_test)
                return cor
        # Return default
        return (0.0, 0.0)
    except Exception as ex:
        # handling errors
        print("error, perform analysis on row", row.name, is_in_pathway)
        print(perform_t_test(components_data,
                             remove_gene_from_pathway(matrix_genes, row.name),
                             background_genes=background_genes_data))
        print("Error", ex)


def single_row_analysis_regression(row, components_data, matrix_genes,
                                   background_genes_data, overall_coef,
                                   overall_intercept, is_in_pathway=False,
                                   components_data_permutated=None):
    """
    Method to perform the logitic regression analysis of a single gene in a single pathway.
    :param row: Row of the gene pathway matrix
    :param components_data: Components data matrix
    :param matrix_genes: Boolean Serie which genes in active or not in the pathway
    :param background_genes_data: Serie with the genes names which must be used
    :param overall_coef: The betas of the overall logistic model analysis between the component values of genes active against genes not active in the pathway
    :param overall_intercept: The intercept of the overall logistic model analysis between the component values of genes active against genes not active in the pathway
    :param is_in_pathway: Boolean if the gene is active in the original pathway description
    :param components_data_permutated: Permutated components data matrix
    """

    # Check if the current gene to analyse is active in the original
    # pathway annotation
    if is_in_pathway:
        #gene is active

        # Remove the gene from the pathway (set the gene to false in the
        # annotation overview) and perform a new logistic regression analysis
        coef, intercept, model = perform_regression(components_data,
                                                    remove_gene_from_pathway(matrix_genes, row.name),
                                        background_genes=background_genes_data)

        # calculate the log2_logodds value
        influence_value = row.T.dot(coef)
        logodds = np.exp(influence_value  + intercept[0])
        logodds_log2 = np.log2(logodds)

        p_val = scipy.special.expit(influence_value + intercept[0])

        # use the permutation matrix to create a specific null distribution for this pathway
        permutated_influence_value = components_data_permutated.dot(coef)
        permutated_logodds = np.exp(permutated_influence_value + intercept[0])
        permutated_logodds_log2 = np.log2(permutated_logodds)
        permutation_mean = permutated_logodds_log2.mean()
        permutation_std = permutated_logodds_log2.std()

        # Calculate the z score
        z_score = (logodds_log2 - permutation_mean) / permutation_std

        return (z_score, p_val)
    else:
        # The gene is not included in the pathway, so the overall logistic
        # regression can be used to calculate the logodds values
        influence_value = row.T.dot(overall_coef)
        logodds = np.exp(influence_value + overall_intercept[0])
        logodds_log2 = np.log2(logodds)

        p_val = scipy.special.expit(influence_value + overall_intercept[0])
        return (logodds_log2, p_val)


def perform_single_pathway(components_data, matrix_genes,
                           background_genes_data=None, analysis_type=None,
                           components_data_permutated=None):
    ###
    ### Function which analysis a single pathway.
    # The input contains the components matrix and a boolean vector of a
    # single pathway. The background genes set can also be given and the
    # analysis type. By default the new regression method will be used to
    # analyse the data but the earlier version of gene network with the
    # t-test approach can also be chosen (analyse_type = t_test).
    # the method returns a series of z-scores for every gene in the pathway.
    # The executed steps are:
    # - Checking which method must be used (regression or t-test)
    ### In case of the Regression analysis:
    # - Perform the overall analysis on all the (background) genes
    # - Use the model to calculate the log2 logodds scores of every gene.
    #   Formula: gene_logg_odd_score = log2(e^(intercept + beta_1 * component_1 + beta_2 * component_2 ... + beta_n * component_n)).
    #   De components scores are extracted from the pathway/components input
    #   Serie which is given as function input (components_data).
    # - Use the model output to calculate the mean and standard error of
    #   the null-distribution by using the permutated matrix.
    # - Use the mean and std of the null-distribution to calculate the
    #   z-scores output of the genes which are not present in the original
    #   pathway annotation (false in the pathway matrix)
    # - For every gene which are present in the original pathway annotation
    #   (true in the pathway matrix), we set that particular gene temporary on
    #   False and fit a new model with the temporary gene annotation information.
    #   We also recalculate the mean and std of the null-distribution with the
    #   new fitted model and calculate for these genes the z-scores.
    # - Combine the z-scores of the genes present and not present in the
    #   original pathway scores (results of the two steps above)
    # - Return the z-scores
    ### In case of the t-test method:
    # - Fit a overall t-test between the component values of genes included in
    #   the pathway annotation and the component values of the genes which are
    #   not annotated in the pathway, which was eventually filtered on
    #   background genes.
    # - The extracted p-value were used to calculate the z-score including the
    #   direction extracted from t-value sign
    # - These overall scores were used as z-score for the genes which are not
    #   present in the pathway annotation matrix (false in the pathway matrix)
    # - For every gene which are present in the original pathway annotation
    #   (true in the pathway matrix), we set that particular gene temporary on
    #   False and perform a new t-test between the component values of genes
    #   included and not included in the pathway.
    # - translate the p-value to a z-score including the direction
    # - Use that z-score as z-scores of that particular gene
    # - The pearson p- and t-value extracted are calculated between the
    #   z-scores and the pathway/components scores to find the gene p-value.
    # - The p-value is translated to a z-score including the direction by using
    #   the sign of the t-value
    # - The z-score of all the genes are returned as a Serie
    ###

    if analysis_type is None or analysis_type == analysis_types["REGRESSION"]:
        # create to dataframes to split up the results of genes which are
        # active in the pathway and the genes which are not active
        # in the pathway
        genes_in_pathway = pd.DataFrame()
        genes_not_in_pathway = pd.DataFrame()

        # Perform the overall regression between the gene/component
        # matrix and the gene pathway annotation
        overall_coefficient, overall_intercept, _ = perform_regression(components_data,
                                              matrix_genes,
                                              background_genes=background_genes_data)
        ###
        ### calculate null distribution based on permutation matrix
        ###

        # apply the model by using the betas and the permutated components matrix
        # Formula: intercept + beta_1 * component_1 + beta_2 * component_2 ... + beta_n * component_n
        overall_permutated_influence_value = overall_intercept[0] + components_data_permutated.dot(overall_coefficient)

        # make the value an exponent of e^(...)
        overall_permutated_logodds = np.exp(overall_permutated_influence_value)

        # Perform log2 transformation
        permutated_logodds_log2 = np.log2(overall_permutated_logodds)

        # Calculate mean and standard deviation
        overall_permutation_mean = permutated_logodds_log2.mean()
        overall_permutation_std = permutated_logodds_log2.std()

        # Calculate the logodds score for all the genes which are not annotated in
        # the original pathway annotation
        genes_not_in_pathway[["r", "p"]] = components_data.loc[
                                           matrix_genes[~matrix_genes].index,
                                           :].apply(
            lambda row: single_row_analysis_regression(row=row,
                                                       components_data=components_data,
                                                       matrix_genes=matrix_genes,
                                                       background_genes_data=background_genes_data,
                                                       overall_coef=overall_coefficient,
                                                       overall_intercept=overall_intercept,
                                                       is_in_pathway=False,
                                                       components_data_permutated=components_data_permutated
                                                       ),
            axis=1, result_type="expand")

        # translate the log-odds score to z-score by using the mean and std
        # received from the permutation strategy
        genes_not_in_pathway["z_score"] = (genes_not_in_pathway["r"] - overall_permutation_mean) / overall_permutation_std

        # calculate the z-scores for the genes which are present in the gene
        # annotation of the pathway
        # the method will fit a new model where the particular gene is set to
        # false and use the permutation matrix to re calculate the mean and std
        # and translate the logodds score to a z-score
        genes_in_pathway[["r", "p"]] = components_data.loc[
                                       matrix_genes[matrix_genes].index,
                                       :].apply(
            lambda row: single_row_analysis_regression(row=row,
                                                       components_data=components_data,
                                                       matrix_genes=matrix_genes,
                                                       background_genes_data=background_genes_data,
                                                       overall_coef=overall_coefficient,
                                                       overall_intercept=overall_intercept,
                                                       is_in_pathway=True,
                                                       components_data_permutated=components_data_permutated),
            axis=1, result_type="expand")
        genes_in_pathway["z_score"] = genes_in_pathway["r"]

        # combine and return the results as Series
        gene_pw_results = genes_in_pathway.append(genes_not_in_pathway)
        return gene_pw_results.loc[:, 'z_score']
    elif analysis_type == analysis_types["T_TEST"]:
        # create to dataframes to split up the results of genes which are
        # active in the pathway and the genes which are not active
        # in the pathway
        genes_in_pathway = pd.DataFrame()
        genes_not_in_pathway = pd.DataFrame()

        # Perform the overall t-test analyses to find the p-values between
        # the genes and the components
        overall_t_test, _, _ = perform_t_test(components_data,
                                              matrix_genes,
                                              background_genes=background_genes_data)

        # Calculate the p and r values of the genes by using Pearson
        # correlation of all the genes which are present in the original
        # pathway annotation
        genes_in_pathway[["r", "p"]] = components_data.loc[
                                       matrix_genes[matrix_genes].index,
                                       :].apply(
            lambda row: single_row_analysis_t_test(row,
                                                   components_data,
                                                   matrix_genes,
                                                   background_genes_data,
                                                   overall_t_test,
                                                   True)
            , axis=1, result_type="expand")

        # Calculate the p and r values of the genes by using Pearson
        # correlation of all the genes which are not present in the original
        # pathway annotation
        genes_not_in_pathway[["r", "p"]] = components_data.loc[
                                           matrix_genes[~matrix_genes].index,
                                           :].apply(
            lambda row: single_row_analysis_t_test(row,
                                                   components_data,
                                                   matrix_genes,
                                                   background_genes_data,
                                                   overall_t_test,
                                                   False),
            axis=1, result_type="expand")


        # combine the dataframes of genes inside and outside
        # the original pathway annotation
        gene_pw_results = genes_in_pathway.append(genes_not_in_pathway)

        # translate the gene-pathway scores to z scores including the direction
        # by using the sign of the r-value
        z_scores_r = scipy.stats.norm.ppf(1 - gene_pw_results.loc[:, 'p'] / 2) * np.sign(gene_pw_results.loc[:, 'r'])
        z_scores_r[np.isnan(z_scores_r)] = 0.0
        return z_scores_r
    else:
        # Raise error if a method name is given which is not implemented
        raise NotImplementedError("Method '{}' is not implemented".format(analysis_type))


def perform_t_test(X, y, background_genes=None, equal_var=False):
    ###
    ### Method to perform the t-test analysis and translate the p-value to
    # a z-score.
    # The method uses two variables, and two optional variables
    # X = the input gene/component matrix.
    # Y = Gene Serie boolean which genes are active in the original pathway
    # definition
    # background_genes = Serie which genes are present in the background
    # gene set (optional)
    # equal_var = boolean, if the t-test must be performed with
    # equal variable (optional, default is false)
    ###

    # Split the genes which are active in the pathway and which are not
    genes_include = X[y]
    genes_not_included= X[~y]

    # Perform the background selection if necessary
    if background_genes is not None:
        bg_intersect_genes = genes_not_included.index.intersection(background_genes)
        genes_not_included = genes_not_included.loc[bg_intersect_genes, :]

    # perform the t-test
    t_value, p_value = scipy.stats.ttest_ind(genes_include,
                                             genes_not_included,
                                             equal_var=equal_var)

    # translate p-value to z-score and add the direction based on the sign
    z_score = scipy.stats.norm.ppf(1 - p_value / 2) * np.sign(t_value)
    return z_score.ravel(), t_value, p_value


def perform_regression(X, y, background_genes=None):
    ###
    ### Method to perform the regression analysis
    # and return the identified intercept and beta values
    # The method uses two variables, and one optional variable
    # X = the input gene/component matrix.
    # Y = Gene Serie boolean which genes are active in the original pathway
    # definition
    # background_genes = Serie which genes are present in the background
    # gene set (optional)
    # equal_var = boolean, if the t-test must be performed with
    # equal variable (optional, default is false)
    ###

    # set the input variables
    X_filtered = X
    y_filtered = y

    # Perform the background selection if necessary
    if background_genes is not None:
        X_filtered = X.loc[background_genes, :]
        y_filtered = y.loc[background_genes]

    #Fit the model
    log_mdl = LogisticRegression(solver="lbfgs", C=1.0,  max_iter=6000, tol=1e-6)
    log_mdl.fit(X_filtered, y_filtered.ravel())

    # return the results
    model_coef = log_mdl.coef_.ravel()
    model_intercept = log_mdl.intercept_
    return model_coef, model_intercept, log_mdl