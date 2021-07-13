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
    create_auc_output_file, remove_gene_from_hpo

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
    def __init__(self, components_data_path, matrix_path,
                 output_dir, analysis_type, minimal_number_of_genes=None,
                 background_genes_path=None, n_cores=None,
                 split_start=None, split_end=None, multi_node_node_id=None,
                 multi_node_num_nodes=None, multi_node_output_dir=None,
                 force=False, output_disable_txt=False,
                 output_disable_gzip=False, output_disable_pickle=False
                 ):
        self.start_time = time.time()
        self.components_data_path = components_data_path
        self.output_dir = output_dir
        self.analysis_type = analysis_type
        self.matrix_path = matrix_path
        self.background_genes_path = background_genes_path

        self.minimal_number_of_genes = -1
        if minimal_number_of_genes is not None and minimal_number_of_genes != "":
            self.minimal_number_of_genes = int(minimal_number_of_genes)
        if self.minimal_number_of_genes < 3:
            self.minimal_number_of_genes = 3
            logging_print("Minimal number of genes must be 3 of higher, so the value is set to 3")

        # for selecting and multiprocessing
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
        self.write_base_loginfo()
        self.read_files()
        self.perform_gene_intersection()
        self.perform_filtering()
        self.perform_multinode_processing()
        self.calculate_gene_pathway_count()
        self.create_permutations_matrixes()
        self.perform_analysis()

        self.calculate_auc_values()
        self.calculate_wilcox_p_value()
        if self.multi_node_num_nodes == None :
            self.calculate_p_value_bonferroni_correction()
        self.write_output_files()
        self.remove_temp_files()

        if self.multi_node_num_nodes != None:
            self.merge_output_files()

        self.write_last_loginfo()

    def write_base_loginfo(self):
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
        timer_print(self.start_time, prefix="## ANALYSIS READY")

    def read_files(self):
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
        print(matrix_data_info)
        logging.info(matrix_data_info)

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
        return self.background_genes_data is not None

    def perform_gene_intersection(self):
        pgi_start_time = time.time()
        intersect_genes = self.components_data.index.intersection(self.matrix_data.index)
        self.components_data = self.components_data.loc[intersect_genes, :]
        self.matrix_data = self.matrix_data.loc[intersect_genes, :]
        if self.has_background_genes():
            self.background_genes_data = pd.Series(list(set(self.background_genes_data).intersection(set(intersect_genes))))

        timer_print(pgi_start_time, prefix="Gene intersection ready")

    def perform_filtering(self):
        if self.minimal_number_of_genes > 0:
            pf_start_time = time.time()
            matrix_selection = self.matrix_data.sum(axis=0) >= self.minimal_number_of_genes
            self.matrix_data  = self.matrix_data.loc[:, matrix_selection[matrix_selection].index]
            logging_print("Minimal gene filtering: {} number of pathways over".format(self.matrix_data.shape[1]))
            timer_print(pf_start_time, prefix="Minimal gene in pathway filtering ready")

    def create_permutations_matrixes(self):
        perm_start_time = time.time()
        permutation_path = "{}_permutation_matrix.pkl".format(self.components_data_path)

        if os.path.isfile(permutation_path) and self.force == False:
            self.components_data_permutated = pd.read_pickle(permutation_path)
            logging_print("Permutation dataframe is loaded from '{}', "
                          "shape: {}".format(permutation_path,
                                             self.components_data_permutated.shape))
        else:
            self.components_data_permutated = self.components_data.set_index(
                np.random.permutation(self.components_data.index))
            self.components_data_permutated.to_pickle(permutation_path)
            logging_print("Permutation table is created, shape: {}".format(
                self.components_data_permutated.shape))
        timer_print(perm_start_time, prefix="permutation is ready")

    def perform_multinode_processing(self):
        # handle new start end stop
        if self.split_end is not None:
            self.matrix_data = self.matrix_data.iloc[:, :self.split_end]

        if self.split_start is not None:
            self.matrix_data = self.matrix_data.iloc[:, self.split_start:]

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
            num_pathways = self.matrix_data.shape[1]
            pathways_per_node = num_pathways // self.multi_node_num_nodes
            if pathways_per_node * self.multi_node_num_nodes < num_pathways:
                pathways_per_node += 1

            start_id = self.multi_node_node_id * pathways_per_node
            end_id = (self.multi_node_node_id + 1) * pathways_per_node
            if end_id > num_pathways:
                end_id = num_pathways

            self.matrix_data = self.matrix_data.iloc[:, start_id:end_id]

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
        self.gene_pathway_count = self.matrix_data.sum(axis=0)

    def perform_analysis(self):
        pa_start_time = time.time()
        total_z_score_results = []

        # check if some temp files where present
        temp_z_score_paths = glob.glob(os.path.join(self.output_dir,
                                                    "temp_results_analysis_z_scores_*.pkl"))
        temp_z_scores = None
        temp_already_processed_pathways = None
        if len(temp_z_score_paths) > 0:
            temp_z_scores = pd.read_pickle(temp_z_score_paths[0])
            logging_print("Temp file '{}' with already processed pathways loaded. size df: {}".format(temp_z_score_paths[0], temp_z_scores.shape))
            temp_already_processed_pathways = list(temp_z_scores.columns.to_numpy())

        pathway_manager = mp.Manager()
        pathway_queue = pathway_manager.Queue()

        retults_manager = mp.Manager()
        results_queue = retults_manager.Queue()

        for index, row in self.matrix_data.iteritems():
            if temp_z_scores is not None:
                if index in temp_already_processed_pathways:
                    total_z_score_results.append(temp_z_scores.loc[:, index])
                else:
                    pathway_queue.put(index)
            else:
                pathway_queue.put(index)

        logging_print("total pathways already done: {}".format(len(total_z_score_results)))

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

        for process in processes:
            process.start()

        total_done = 0
        last_save_time = time.time()
        save_time_in_minutes = 2
        z_score_file_path = None
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
                        temp_dataframe_z_scores.to_pickle(new_z_score_file_path)

                        if z_score_file_path is not None and os.path.isfile(z_score_file_path):
                            os.remove(z_score_file_path)
                        z_score_file_path = new_z_score_file_path

            except queue.Empty:
                time.sleep(1)
                continue

        for process in processes:
            process.join()

        print("all processes ready")
        self.pathway_gene_scores = pd.DataFrame(total_z_score_results).T

        pathway_gene_scores_temp_file_path = os.path.join(self.output_dir,
                                             "temp_pathway_gene_scores_temp_{}.pkl".format(
                                                 datetime.now()))
        self.pathway_gene_scores.to_pickle(pathway_gene_scores_temp_file_path)
        timer_print(pa_start_time, prefix="Performing analysis ready")

    def calculate_auc_values(self):
        cav_st_stime = time.time()
        auc_values = {}
        for pathway_id in self.pathway_gene_scores.columns:
            if pathway_id in list(self.matrix_data.columns):
                pathway_gene_set = self.matrix_data.loc[:, pathway_id]
                pathway_gene_scores = self.pathway_gene_scores.loc[self.matrix_data.index, pathway_id]
                try:
                    log_reg_fpr, log_reg_tpr, _ = roc_curve(pathway_gene_set * 1,
                                                            pathway_gene_scores)
                    auc_values[pathway_id] = auc(log_reg_fpr, log_reg_tpr)
                except ValueError:
                    logging_print("Value error in auc calculation of pathway: {}".format(pathway_id))
            else:
                logging_print("HPO term: {} not in matrix".format(pathway_id))

        auc_calc_values = pd.Series(auc_values)
        auc_calc_values[auc_calc_values > 1.0] = 1.0
        self.auc_values = auc_calc_values
        timer_print(cav_st_stime, prefix="AUC value calculation ready")

    def calculate_wilcox_p_value(self):
        cwpv_st_time = time.time()
        wal_p_values = {}
        for pathway_id in self.auc_values.index:
            if pathway_id in list(self.matrix_data.columns):
                selected_pathway = self.pathway_gene_scores.loc[self.matrix_data.index, pathway_id]
                include_in_pathway = selected_pathway[self.matrix_data.loc[:, pathway_id]]
                include_not_in_pathway = selected_pathway[~self.matrix_data.loc[:, pathway_id]]
                try:
                    _, p_value = scipy.stats.mannwhitneyu(include_in_pathway,
                                                       include_not_in_pathway,
                                                       use_continuity=True,
                                                       alternative="two-sided")
                    wal_p_values[pathway_id] = p_value
                except ValueError:
                    logging_print("Value error in p value calculation, pvalue 0.0 is set to pathway: {}".format(pathway_id))
                    wal_p_values[pathway_id] = 0.0

        self.p_values = pd.Series(wal_p_values)
        timer_print(cwpv_st_time, prefix="Wilcox p value calculation ready")

    def calculate_p_value_bonferroni_correction(self, alpha=0.05):
        cpvbc_st_time = time.time()
        self.bonf_p_values = self.p_values * self.p_values.shape[0]
        timer_print(cpvbc_st_time,
                         prefix="Bonferroni p value correction ready")

    def write_output_files(self):
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
        for file_path in glob.glob(os.path.join(self.output_dir,
                                                "temp_results_analysis_z_scores_*.pkl")):
            os.remove(file_path)

    def merge_output_files(self):
        mof_start_time = time.time()
        if self.multi_node_output_dir:
            node_output_auc_pred_files = glob.glob(os.path.join(
                self.multi_node_output_dir, "*", "predictions_auc.pkl"))
            if len(node_output_auc_pred_files) == self.multi_node_num_nodes:
                logging_print("Merge output files")

                # merge AUC prediction files
                comp_df_auc_list = []
                for node_output_auc_file_path in node_output_auc_pred_files:
                    file_auc_df = pd.read_pickle(node_output_auc_file_path)
                    comp_df_auc_list.append(file_auc_df)

                overall_auc_df = pd.concat(comp_df_auc_list)

                bonf_values = overall_auc_df["pValue"] * overall_auc_df.shape[0]
                bonf_values[bonf_values > 1] = 1
                bonf_values[bonf_values < 0] = 0
                overall_auc_df["bonferroni"] = bonf_values

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

                if self.output_disable_pickle == False:
                    overall_auc_df.to_pickle(os.path.join(
                        self.multi_node_output_dir, "predictions_auc_bonf.pkl"))

                # merge gene pathway files
                node_output_gene_pathway_pred_files = glob.glob(
                    os.path.join(self.multi_node_output_dir, "*",
                                 "gene_pathway_scores.pkl"))
                comp_df_gene_pathway_list = []
                for node_output_auc_file_path in node_output_gene_pathway_pred_files:
                    file_gene_pathway_df = pd.read_pickle(node_output_auc_file_path)
                    comp_df_gene_pathway_list.append(file_gene_pathway_df)

                overall_gene_pathway_df = pd.concat(comp_df_gene_pathway_list,
                                                    axis=1)

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

                if self.output_disable_pickle == False:
                    overall_gene_pathway_df.to_pickle(
                        os.path.join(self.multi_node_output_dir,
                                     "gene_pathway_scores.pkl"))
        timer_print(mof_start_time, prefix="Writing merged outputfile ready")


def single_pathway_worker(components_data, matrix_data, pathway_queue,
                          results_queue, background_genes_data=None,
                          analysis_type=None, max_end_time=-1,
                          components_data_permutated=None):
    try:
        while True:
            try:
                if pathway_queue.empty():
                    break
                if 0 < max_end_time <= time.time():
                    break

                single_pathway_id = pathway_queue.get(True, timeout=1)
                results = perform_single_pathway(components_data,
                                                 matrix_data.loc[:, single_pathway_id],
                                                 background_genes_data,
                                                 analysis_type,
                                                 components_data_permutated)
                results.name = single_pathway_id
                results_queue.put(results)
            except queue.Empty:
                time.sleep(1)
                continue

        # worker is done
        results_queue.put("DONE")
        time.sleep(1)
        return

    except Exception as ex:
        print("ERROR in worker: {}".format(ex))
        print("exception done")
        results_queue.put("DONE")
        time.sleep(1)
        logging.exception("ERROR in worker")


def single_row_analysis_t_test(row, components_data,
                               matrix_genes, background_genes_data,
                               overall_test, is_in_pathway=False):
    try:
        if is_in_pathway:
            new_gene_set = remove_gene_from_hpo(matrix_genes, row.name)
            intersect_components_data_keys = components_data.index.intersection(new_gene_set.index)

            z_scores, _, _ = perform_t_test(components_data.loc[intersect_components_data_keys, :],
                                            new_gene_set,
                                            background_genes=background_genes_data)
            if z_scores.sum() > 0:
                cor = scipy.stats.pearsonr(row, z_scores)
                return cor
        else:
            if overall_test.sum() > 0:
                cor = scipy.stats.pearsonr(row, overall_test)
                return cor
        return (0.0, 0.0)
    except Exception as ex:
        print("error, perform analysis on row", row.name, is_in_pathway)
        print(perform_t_test(components_data,
                             remove_gene_from_hpo(matrix_genes, row.name),
                             background_genes=background_genes_data))
        print("Error", ex)


def single_row_analysis_regression(row, components_data, matrix_genes, background_genes_data, overall_coef, overall_intercept, is_in_pathway=False, components_data_permutated=None):
    if is_in_pathway:
        coef, intercept, model = perform_regression(components_data,
                                        remove_gene_from_hpo(matrix_genes, row.name),
                                        background_genes=background_genes_data)

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
        z_score = (logodds_log2 - permutation_mean) / permutation_std

        return (z_score, p_val)
    else:
        influence_value = row.T.dot(overall_coef)
        logodds = np.exp(influence_value + overall_intercept[0])
        logodds_log2 = np.log2(logodds)

        p_val = scipy.special.expit(influence_value + overall_intercept[0])
        return (logodds_log2, p_val)


def perform_single_pathway(components_data, matrix_genes, background_genes_data=None, analysis_type=None, components_data_permutated=None):
    if analysis_type is None or analysis_type == analysis_types["REGRESSION"]:
        genes_in_pathway = pd.DataFrame()
        genes_not_in_pathway = pd.DataFrame()

        # Perform overall regression
        overall_coefficient, overall_intercept, _ = perform_regression(components_data,
                                              matrix_genes,
                                              background_genes=background_genes_data)

        # calculate null distribution based on permutation matrix
        overall_permutated_influence_value = components_data_permutated.dot(overall_coefficient)
        overall_permutated_logodds = np.exp(overall_permutated_influence_value + overall_intercept[0])
        permutated_logodds_log2 = np.log2(overall_permutated_logodds)
        overall_permutation_mean = permutated_logodds_log2.mean()
        overall_permutation_std = permutated_logodds_log2.std()
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
        genes_not_in_pathway["z_score"] = (genes_not_in_pathway["r"] - overall_permutation_mean) / overall_permutation_std
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
        gene_pw_results = genes_in_pathway.append(genes_not_in_pathway)
        return gene_pw_results.loc[:, 'z_score']
    elif  analysis_type == analysis_types["T_TEST"]:
        genes_in_pathway = pd.DataFrame()
        genes_not_in_pathway = pd.DataFrame()
        overall_t_test, _, _ = perform_t_test(components_data,
                                              matrix_genes,
                                              background_genes=background_genes_data)

        genes_in_pathway[["r", "p"]] = components_data.loc[matrix_genes[matrix_genes].index, :].apply(lambda row: single_row_analysis_t_test(row, components_data, matrix_genes, background_genes_data, overall_t_test, True), axis=1, result_type="expand")
        genes_not_in_pathway[["r", "p"]] = components_data.loc[matrix_genes[~matrix_genes].index, :].apply(lambda row: single_row_analysis_t_test(row, components_data, matrix_genes, background_genes_data, overall_t_test, False), axis=1, result_type="expand")
        gene_pw_results = genes_in_pathway.append(genes_not_in_pathway)
        z_scores_r = scipy.stats.norm.ppf(1 - gene_pw_results.loc[:, 'p'] / 2) * np.sign(gene_pw_results.loc[:, 'r'])
        z_scores_r[np.isnan(z_scores_r)] = 0.0
        return z_scores_r
    else:
        raise NotImplementedError("Method '{}' is not implemented".format(analysis_type))


def perform_t_test(X, y, background_genes=None, equal_var=False):
    genes_include = X[y]
    genes_not_included= X[~y]
    if background_genes is not None:
        bg_intersect_genes = genes_not_included.index.intersection(background_genes)
        genes_not_included = genes_not_included.loc[bg_intersect_genes, :]

    t_value, p_value = scipy.stats.ttest_ind(genes_include,
                                             genes_not_included,
                                             equal_var=equal_var)

    z_score = scipy.stats.norm.ppf(1 - p_value / 2) * np.sign(t_value)
    return z_score.ravel(), t_value, p_value


def perform_regression(X, y, background_genes=None):
    X_filtered = X
    y_filtered = y
    if background_genes is not None:
        X_filtered = X.loc[background_genes, :]
        y_filtered = y.loc[background_genes]

    log_mdl = LogisticRegression(solver="lbfgs", C=1.0,  max_iter=6000, tol=1e-6)
    log_mdl.fit(X_filtered, y_filtered.ravel())
    model_coef = log_mdl.coef_.ravel()
    model_intercept = log_mdl.intercept_
    return model_coef, model_intercept, log_mdl