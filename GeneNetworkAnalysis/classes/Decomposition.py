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
    ###
    ### Class to perform the decomposition part of the gene network analysis
    ###
    def __init__(self, input_file_path, output_dir, analysis_type,
                 pca_type=None, over_samples=False, n_components=None,
                 test_run=False, fastICA_max_iter=None,
                 fastICA_stable_iterations=None,
                 fastica_stable_safe_intermediates=False, n_rows=None,
                 perform_log2=False, preprocessing_center_scale=False,
                 force=False, output_disable_txt=False,
                 output_disable_gzip=False, output_disable_pickle=False):
        """
        The input of this class is comparable with the input variables of the script
        :param input_file_path: The path of the input file. The input file contains the gene expression values where the rows correspond to the samples and the columns to the genes. The input file can be given as pandas pickled dataframe or as tab separated txt file. The txt input file is automatically decompressed if the file paths ends withs ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’.
        :param output_dir: Path to the output directory. The output files will be created in this directory. The directory will automatically be created if not exists.
        :param analysis_type: Set the analysis type which will be used to perform the decomposition. Can be set to ‘PCA’, ‘FASTICA’, ‘FASTICA_STABLE’ or ‘FASTICA_STABLE_FROM_TEMP’ PCA: Use principal component analysis as decomposition method. <br><br> PCA is implemented by using SVD algorithms Different algorithms are implemented and can be selected by using the ‘pca_type’ parameter. <br><br> FASTICA: Use the fastICA implementation to perform an independent component analysis as decomposition method. The whiting step is done by the implemented SVD algorithms and can be selected by using the ‘pca_type’ parameter. <br><br> FASTICA_STABLE: Use multiple different fastICA runs and calculate the average gene eigenvector loadings of the independent components to get sable independent components. The FastICA algorithms use a random initialization which can result in slightly different outcomes per individual run. This analysis can take a long time. <br><br> FASTICA_STABLE_FROM_TEMP: Method to finish up an uncomplete “FASTICA_STABLE” run. The executed runs will be combined to the final independent components. The method will not execute more FastCIA runs even if the number of set iterations is not reached (the option: fastica_stable_iterations)
        :param pca_type: Select the SVD type. The the analysis type is set to FastICA, the selected PCA type will be used to performing the whiting step.
        :param over_samples: Transpose the input matrix before the analysis. The columns will be processed as rows and the rows as columns. The output matrixes are not changed.
        :param n_components: The number of components which will be returned in the export files. If the number is not set, the maximal number of components will be returned
        :param test_run: Perform a test run with the first 150 samples and 100 genes
        :param fastICA_max_iter: Maximal number of FastICA iterations. Can only be used if analysis is set to FastICA or FastICA-stable. default is 2000
        :param fastICA_stable_iterations: Number of FastICA runs used to calculate the final FastICA stable results. Can only be set if analysis is set to FastICA-stable. Default is 10
        :param fastica_stable_safe_intermediates: Safe the intermediate FastICA runs in the fastICA-stable run. Can only be set if analysis is set to FastICA-stable.
        :param n_rows: Number of rows in the input dataset to process.
        :param perform_log2: Perform log2 transformation of the input data before the analysis
        :param preprocessing_center_scale: Perform center scale transformation of the input data before the analysis
        :param force: Do not use the cached versions of the input data but the orginal input files.
        :param output_disable_txt: Do not export the output as text files.
        :param output_disable_gzip: Do not gzip the data output text files
        :param output_disable_pickle: Do not export the output as pickle Pandas files.
        """

        # Set the variables
        self.input_file_path = input_file_path
        self.output_dir = output_dir
        self.analysis_type = analysis_type
        self.n_components = None
        if n_components is not None and n_components != "":
            self.n_components = int(n_components)

        # set default PCA type
        self.pca_type = "auto"
        if pca_type is not None:
            self.pca_type = pca_type

        self.over_samples = over_samples
        self.test_run = test_run

        # set default number of FastICA iterations
        self.fastICA_max_iter = 2000
        if fastICA_max_iter != None and fastICA_max_iter != "":
            self.fastICA_max_iter = fastICA_max_iter

        # set default number of FastICA stable iterations
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

        # set output params
        self.output_disable_txt = output_disable_txt
        self.output_disable_gzip = output_disable_gzip
        self.output_disable_pickle = output_disable_pickle

        # create results output variables
        self.projected_data = None
        self.components = None
        self.method_stats = {}

    def run(self, create_output_files=True,
            write_end_of_log_file=True ,
            output_files_prefix=None):
        ###
        ### method to run the complete analysis
        ###

        # Read input data if input data is not already present
        if self.input_data is None:
            # create logfile and write the basic analysis log information
            self.write_base_loginfo()

            # read all the input files
            self.read_file()

            # perform the preprocessing steps
            self.perform_data_preprocessing()

            # Create the analysis opject based on the input parameters.
            # The analysis objects uses the same interface, so after creating
            # the analysis object based on the input parameters the objects can
            # be used in the same way without knowing which type of analysis
            # will be performed. This makes it possible to add
            # more decomposition algorithms in the future
            self.set_analysis_object()

        # Perform the analysis
        self.perform_analysis()

        # create the output files
        if create_output_files:
            self.write_output_files(file_prefix=output_files_prefix)

        # Close the logfile
        if write_end_of_log_file:
            self.write_last_loginfo()
            self.write_timeit_export_file()

    def write_base_loginfo(self):
        #
        # Write the basis log information
        #
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
        #
        # Add the end information to the logfile
        #
        timer_print(self.start_time, prefix="## ANALYSIS READY",
                    time_overview_log=self.process_time_overview)

    def read_file(self):
        #
        # Method to read all the input files
        #

        # save the start time
        rf_start_time = time.time()

        # Read the input file if present
        if os.path.isfile(self.input_file_path):
            if self.test_run:
                # Read a test dataset, so only the first 150 columns and rows
                self.input_data = pd.read_csv(self.input_file_path,
                                              sep="\t", index_col=0,
                                              nrows=150)
                self.input_data = self.input_data.iloc[:150, :100]
            else:
                # Check if only a part of the rows must be loaded instead
                # of the complete dataset
                n_rows = None
                if self.n_rows != None and self.n_rows != '':
                    n_rows = int(self.n_rows)

                # read the input file, This can be a (cahsed) Pandas pickle
                # file or a tab seperated text file which can be compressed.
                # If Force is set to true, the method will not read cashed
                # pickle files created from the original txt matrix if these
                # are present. By default the cashed version will be
                # loaded (which contain the suffix _cashed.pickle) if this file
                # is present in the same directory as the input matrix

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

        # if the number of components is not set, we will set this to
        # the smallest direction of the input matrix
        if self.n_components is None:
            self.n_components = np.min(self.input_data.shape)

        # log some basic info
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
        # get the matrix used for the analysis.
        # The matrix is transposed if the parameter "over samples"
        # is set to True
        return self.input_data if self.over_samples else self.input_data.T

    def perform_data_preprocessing(self):
        #
        # Perform the preprocessing steps
        #
        pp_start_time = time.time()

        # perform log2 transformation over the input data
        if self.perform_log2:
            self.input_data = np.log2(self.input_data,
                                      where=(self.input_data != 0.0))
            logging_print("Log 2 transformation is performed on the input data")

        # perform center scaling over the input data
        if self.pre_processing_center_scale:
            scaler = StandardScaler()
            if self.over_samples:
                scaled_np_data = scaler.fit_transform(self.input_data)
            else:
                # Transform the dataset before and after scaling
                scaled_np_data = scaler.fit_transform(self.input_data.T).T

            # The method return a Numpy dataframe. Create a new Pandas DataFrame
            self.input_data = pd.DataFrame(scaled_np_data,
                                           index=self.input_data.index,
                                           columns=self.input_data.columns)
            logging_print("Centering and scaling is performed on the input data")

        timer_print(pp_start_time,
                    prefix="Data pre-processing ready",
                    time_overview_log=self.process_time_overview)

    def set_analysis_object(self):
        ##
        ## This mehtod creates the analysis object (wrapper objects) based on
        # the input parameters. All the analysis objects implements the same
        # interface, so after creating the object, the analysis can be
        # performed in the same way for every decomposition type which makes it
        # also easy to add new decompositions in the future.
        ##

        # Create FastICA decomposition object
        if self.analysis_type == decomposition_types["FASTICA"]:
            self.analysis_object = FastICA_wrapper(
                svd_type=self.pca_type,
                n_components=self.n_components,
                max_iter=self.fastICA_max_iter
            )

        # Create PCA decomposition object
        elif self.analysis_type == decomposition_types["PCA"]:
            self.analysis_object = SVD_wrapper(svd_type=self.pca_type,
                                               output_dir=self.output_dir,
                                  n_components=self.n_components)
        # Create stabilised version of FastICA decomposition object
        # This method performs multiple normal FastICA decompositions and
        # calculates the centriods between the eigenvectors of
        # each run. The Centriod object uses a FastICA object as input which
        # will be used in the individual FastICA runs
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

        # Same as FastICA stable but the results will be loaded from an ealier
        # run and than the centriods will be calculated. So this analysis will
        # not run FastICA decompositions but will only finish an earlier
        # (crashed) analysis
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
        #
        # Method to perform an analysis based on the earlier created analysis object
        #
        pa_start_time = time.time()

        # All analysis classes implements an fit method which uses
        # the given input matrix to perform the decomposition
        self.analysis_object.fit(self.get_analysis_input_df())

        # print the used time to the logging file
        timer_print(pa_start_time,
                    prefix="Performing analysis ready",
                    time_overview_log=self.process_time_overview)


    def write_output_files(self, file_prefix=None):
        #
        # Method to create the output files
        #

        wof_start_time = time.time()

        # create the output eigenvector and pc scores matrix,
        # which will contains always the same information.
        # Normally the eigenvectors contains the components
        # (genes / components matrix) and the pc scores contains the rotated data (sample
        # / components) matrix. if the over_samples parameter is set,
        # the eigenvectors file contains the rotated data
        # (genes / components matrix) and the pc scores contains the
        # components (sample / components) matrix.

        if self.over_samples is False:
            eigenvectors_df = self.analysis_object.components.T
            pc_scores_df = self.analysis_object.projected_data
        else:
            eigenvectors_df = self.analysis_object.projected_data
            pc_scores_df = self.analysis_object.components.T

        # add the date information to the index name (first column of the
        # output files)
        eigenvectors_df.index.name = datetime.now().strftime('%d/%m/%Y')
        pc_scores_df.index.name = datetime.now().strftime('%d/%m/%Y')

        # Create the file names
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

        # wirte the export time to the log file
        timer_print(wof_start_time,
                    prefix="Writing output files ready",
                    time_overview_log=self.process_time_overview)

        # if fastICA Stable analysis is used and the option to safe the
        # intermediate steps are set, the code below will export these data
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


# The different decomposition options defined as definitions
decomposition_types = {
    "FASTICA": "fastica",
    "PCA": "pca",
    "FASTICA_STABLE": "fastica-stable",
    "FASTICA_STABLE_FROM_TEMP": "fastica-stable-temp"
}