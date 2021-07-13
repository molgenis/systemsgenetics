"""
File:         Decomposition.py
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

# enable faulthandler to get a trace of Segmentation fault errors
import faulthandler
faulthandler.enable()

import os
import time
import numpy as np
import pandas as pd
import json
from datetime import datetime
from sklearn.preprocessing import StandardScaler
from classes.FastICA_wrapper import FastICA_wrapper
from classes.Centroid_decomposition_wrapper import Centroid_decomposition_wrapper
from classes.SVD_wrapper import SVD_wrapper
from utilities import read_pd_df, stats_dict_to_string, \
    logging_print, timer_print, create_output_dir_if_not_exists


class Decomposition:
    def __init__(self, input_file_path, output_dir, analysis_type,
                 pca_type=None, over_samples=False, n_components=None,
                 test_run=False, fastICA_max_iter=None,
                 fastICA_stable_iterations=None,
                 fastica_stable_safe_intermediates=False, n_rows=None,
                 perform_log2=False, preprocessing_center_scale=False,
                 force=False, output_disable_txt=False,
                 output_disable_gzip=False, output_disable_pickle=False):

        self.input_file_path = input_file_path
        self.output_dir = output_dir
        self.analysis_type = analysis_type
        self.n_components = None
        if n_components is not None and n_components != "":
            self.n_components = int(n_components)

        self.pca_type = "auto"
        if pca_type is not None:
            self.pca_type = pca_type

        self.over_samples = over_samples
        self.test_run = test_run
        self.fastICA_max_iter = 2000
        if fastICA_max_iter != None and fastICA_max_iter != "":
            self.fastICA_max_iter = fastICA_max_iter

        self.fastICA_stable_iterations = 10
        if fastICA_stable_iterations != None and fastICA_stable_iterations != "":
            self.fastICA_stable_iterations = fastICA_stable_iterations

        self.fastica_stable_safe_intermediates = fastica_stable_safe_intermediates
        self.n_rows = n_rows
        self.force = force
        self.perform_log2 = perform_log2
        self.pre_processing_center_scale = preprocessing_center_scale
        self.start_time = time.time()
        self.process_time_overview = {}

        self.input_data = None
        self.analysis_object = None

        # output
        self.output_disable_txt = output_disable_txt
        self.output_disable_gzip = output_disable_gzip
        self.output_disable_pickle = output_disable_pickle

        # results
        self.projected_data = None
        self.components = None
        self.method_stats = {}

    def run(self, create_output_files=True,
            write_end_of_log_file=True ,
            output_files_prefix=None):
        if self.input_data is None:
            self.write_base_loginfo()
            self.read_file()
            self.perform_data_preprocessing()
            self.set_analysis_object()
        self.perform_analysis()

        if create_output_files:
            self.write_output_files(file_prefix=output_files_prefix)
        if write_end_of_log_file:
            self.write_last_loginfo()
            self.write_timeit_export_file()

    def write_base_loginfo(self):
        logging_print(stats_dict_to_string({
            "## START DECOMPOSITION ##": "",
            "DATE": datetime.now(),
            "Input file": self.input_file_path,
            "Output dir": self.output_dir,
            "Analysis type": self.analysis_type,
            "Over samples": self.over_samples,
            "Test run": self.test_run,
            "FastICA max iter": self.fastICA_max_iter,
            "Number of components": self.n_components,
            "Number of rows": self.n_rows,
            "Perform log2 transformation": self.perform_log2,
            "Perform centering and scaling": self.pre_processing_center_scale,
            "Force": self.force
        }))

    def write_last_loginfo(self):
        timer_print(self.start_time, prefix="## ANALYSIS READY",
                    time_overview_log=self.process_time_overview)

    def read_file(self):
        rf_start_time = time.time()

        if os.path.isfile(self.input_file_path):
            if self.test_run:
                self.input_data = pd.read_csv(self.input_file_path,
                                              sep="\t", index_col=0,
                                              nrows=150)
                self.input_data = self.input_data.iloc[:150, :100]
            else:
                n_rows = None
                if self.n_rows != None and self.n_rows != '':
                    n_rows = int(self.n_rows)
                self.input_data = read_pd_df(
                    self.input_file_path,
                    {
                        "sep": "\t",
                        "index_col": 0,
                        "nrows": n_rows
                    },
                    force=self.force)
        else:
            raise FileNotFoundError("Cannot find input file: {}".format(
                self.input_file_path
            ))

        if self.n_components is None:
            self.n_components = np.min(self.input_data.shape)

        logging_print(stats_dict_to_string({
            "Input dataframe n_row": self.input_data.shape[0],
            "Input dataframe n_col": self.input_data.shape[1],
            "first column headers": self.input_data.columns.values[:5],
            "first row index": self.input_data.index.values[:5]
        }))
        timer_print(rf_start_time,
                    prefix="Reading input file ready",
                    time_overview_log=self.process_time_overview)

    def get_analysis_input_df(self):
        return self.input_data if self.over_samples else self.input_data.T

    def perform_data_preprocessing(self):
        pp_start_time = time.time()
        if self.perform_log2:
            self.input_data = np.log2(self.input_data,
                                      where=(self.input_data != 0.0))
            logging_print("Log 2 transformation is performed on the input data")

        if self.pre_processing_center_scale:
            scaler = StandardScaler()
            if self.over_samples:
                scaled_np_data = scaler.fit_transform(self.input_data)
            else:
                scaled_np_data = scaler.fit_transform(self.input_data.T).T

            self.input_data = pd.DataFrame(scaled_np_data,
                                           index=self.input_data.index,
                                           columns=self.input_data.columns)
            logging_print("Centering and scaling is performed on the input data")

        timer_print(pp_start_time,
                    prefix="Data pre-processing ready",
                    time_overview_log=self.process_time_overview)

    def set_analysis_object(self):
        if self.analysis_type == decomposition_types["FASTICA"]:
            self.analysis_object = FastICA_wrapper(
                svd_type=self.pca_type,
                n_components=self.n_components,
                max_iter=self.fastICA_max_iter
            )

        elif self.analysis_type == decomposition_types["PCA"]:
            self.analysis_object = SVD_wrapper(svd_type=self.pca_type,
                                               output_dir=self.output_dir,
                                  n_components=self.n_components)

        elif self.analysis_type == decomposition_types["FASTICA_STABLE"]:
            self.analysis_object = Centroid_decomposition_wrapper(
                FastICA_wrapper(
                    svd_type=self.pca_type,
                    n_components=self.n_components,
                    max_iter=self.fastICA_max_iter
                ),
                n_iterations = self.fastICA_stable_iterations,
                output_path=self.output_dir
            )
        elif self.analysis_type == decomposition_types["FASTICA_STABLE_FROM_TEMP"]:
            self.analysis_object = Centroid_decomposition_wrapper(
                FastICA_wrapper(
                    svd_type=self.pca_type,
                    n_components=self.n_components,
                    max_iter=self.fastICA_max_iter
                ),
                n_iterations=self.fastICA_stable_iterations,
                output_path=self.output_dir,
                from_temp_runs=True
            )

    def perform_analysis(self):
        pa_start_time = time.time()
        self.analysis_object.fit(self.get_analysis_input_df())
        timer_print(pa_start_time,
                    prefix="Performing analysis ready",
                    time_overview_log=self.process_time_overview)

    def write_output_files(self, file_prefix=None):
        wof_start_time = time.time()
        if self.over_samples is False:
            eigenvectors_df = self.analysis_object.components.T
            pc_scores_df = self.analysis_object.projected_data
        else:
            eigenvectors_df = self.analysis_object.projected_data
            pc_scores_df = self.analysis_object.components.T
        eigenvectors_df.index.name = datetime.now().strftime('%d/%m/%Y')
        pc_scores_df.index.name = datetime.now().strftime('%d/%m/%Y')

        eigenvectors_file_name = "eigenvectors"
        pc_scores_file_name = "pc-scores"
        if file_prefix is not None:
            eigenvectors_file_name = "{}_{}".format(file_prefix,
                                                    eigenvectors_file_name)
            pc_scores_file_name = "{}_{}".format(file_prefix,
                                                 pc_scores_file_name)

        # export eigenvectors (components)
        if self.output_disable_txt == False:
            if self.output_disable_gzip:
                eigenvectors_df.to_csv(os.path.join(self.output_dir,
                                                    eigenvectors_file_name + ".txt"),
                                       sep='\t')
            else:
                eigenvectors_df.to_csv(os.path.join(self.output_dir,
                                                    eigenvectors_file_name + ".txt.gzip"),
                                       sep='\t')
        if self.output_disable_pickle == False:
            eigenvectors_df.to_pickle(
                os.path.join(self.output_dir, eigenvectors_file_name + ".pkl"))

        # export PC scores
        if self.output_disable_txt == False:
            if self.output_disable_gzip:
                pc_scores_df.to_csv(os.path.join(self.output_dir,
                                                 pc_scores_file_name + ".txt"),
                                    sep='\t')
            else:
                pc_scores_df.to_csv(os.path.join(self.output_dir,
                                                 pc_scores_file_name + ".txt.gzip"),
                                    sep='\t')

        if self.output_disable_pickle == False:
            pc_scores_df.to_pickle(os.path.join(self.output_dir,
                                                pc_scores_file_name + ".pkl"))

        timer_print(wof_start_time,
                    prefix="Writing output files ready",
                    time_overview_log=self.process_time_overview)

        if self.fastica_stable_safe_intermediates and self.analysis_type == decomposition_types["FASTICA_STABLE"]:
            individual_run_dir = os.path.join(self.output_dir,
                                              "FastICA_individual_component_runs")
            create_output_dir_if_not_exists(individual_run_dir)
            for index, run_df in enumerate(self.analysis_object.overall_components):
                ind_df_path = os.path.join(individual_run_dir,
                                           "fastica_components_run_{}.pkl".format(index + 1))
                run_df.to_pickle(ind_df_path)


    def write_timeit_export_file(self):
        # export time overview
        json_file = open(os.path.join(self.output_dir,
                               "processing_time_overview.json"), 'w')
        json.dump(self.process_time_overview, json_file)
        json_file.close()

decomposition_types = {
    "FASTICA": "fastica",
    "PCA": "pca",
    "FASTICA_STABLE": "fastica-stable",
    "FASTICA_STABLE_FROM_TEMP": "fastica-stable-temp"
}