"""
File:         Centroid_decomposition_wrapper.py
Created:      2020-05-05
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

import os
import glob
import time
import numpy as np
import pandas as pd
from utilities import logging_print, timer_print


class Centroid_decomposition_wrapper:
    """
    Class to perform multiple (FastICA) decompositions and calculate the
    average components of all these individual runs. The class is used to
    create a stable decomposition.
    """
    def __init__(self, analysis_object = None, n_iterations=10, output_path=None, from_temp_runs=False):
        """
        The
        :param analysis_object: analysis object which must be fit several times and the resulted components will be used to calculte the average components
        :param n_iterations: int, number of decomposition iterations
        :param output_path: analysis output directory
        :param from_temp_runs: bool, true: do not perform the iterations, but use the earlier excecuted interations to calculate the average components
        """
        self.analysis_object = analysis_object
        self.n_iterations = n_iterations
        self.n_components = analysis_object.n_components
        self.overall_components = []
        self.projected_data = None
        self.components = None
        self.output_path = output_path
        self.from_temp_runs = from_temp_runs

    def fit(self, data):
        # Method to fit the stabilised decomposition

        # List to save the components of the individual decompositions
        self.overall_components = []
        if self.from_temp_runs:
            # Create the average components earlier executed decompositions
            # find the paths to the earlier components
            path = os.path.join(self.output_path,
                                "components_run_*.pickle")
            files_paths = glob.glob(path)

            # Read the earlier files
            for file_path in files_paths:
                comp_run_data = pd.read_pickle(file_path)
                self.overall_components.append(comp_run_data)
            # Calculate the average components
            self.fit_centroids(data)
        else:
            # Run the decompositions
            if self.analysis_object is not None:
                for index in range(self.n_iterations):
                    #  fit the decomposition
                    self.analysis_object.fit(data)
                    # save the components
                    self.overall_components.append(self.analysis_object.components)
                    # export the individual run
                    if self.output_path is not None:
                        path = os.path.join(self.output_path,
                                            "components_run_{}.pickle".format(index+1))
                        self.analysis_object.components.to_pickle(path)
                # Calculate the average components
                self.fit_centroids(data)

    def fit_centroids(self, data):
        # Method to fit the centriods, e.q. calculate the average components
        # from the different runs

        # calculate the average components
        self.components = self.calculate_centriods()
        # Project the on the new components
        proj_data_np =  np.dot(self.components, data.T)
        proj_data_df = pd.DataFrame(proj_data_np,
                                   index=self.overall_components[0].index,
                                   columns=data.index)
        self.projected_data = proj_data_df

    def change_components_number(self, n_components):
        # change the number of components
        self.analysis_object.change_components_number(n_components)
        self.n_components = self.analysis_object.n_components

    def fit_transform(self, data):
        #  Method to fit the decompositions, calculate the average components
        #  and return the projected data
        self.fit(data)
        return self.projected_data

    def calculate_centriods(self):
        # Method to calculate the centriods, e.q. the average components from
        # the components extracted from different individual runs
        calc_cent_start_time = time.time()
        logging_print("Start calculation of the centriods")
        # start position, extracted from the first run
        initial_df = self.overall_components[0]
        gene_names = initial_df.columns
        n_components = initial_df.shape[0]
        initial_np = initial_df.to_numpy()
        signs = np.sign(initial_np)
        initial_np = np.abs(initial_np)
        processed_dfs = [initial_np * signs]
        # loop through all the component runs
        for ica_df in self.overall_components[1:]:
            # Sort the components based on the highest correlation between the
            # components and the component from the first run. this resulted in
            # the same order of the fastica components over all runs.
            # By default, the components of a fastICA are random ordered
            ica_df = ica_df.loc[:, gene_names]
            ica_df = ica_df.iloc[:n_components, :]
            # calculate the correlation
            corr_table = np.abs(
                np.corrcoef(initial_np,
                            np.abs(ica_df.to_numpy())))
            corr_table = corr_table[n_components:, :n_components]

            # find the order based on the maximal correlation
            max_indexes_columns = np.argmax(corr_table, axis=0)
            # Order the components
            ica_df = ica_df.iloc[max_indexes_columns, :]
            # Save the components in the right order
            processed_dfs.append(ica_df.to_numpy())

        # Calculate the centriods e.q. average components
        centriods_array = np.c_[np.abs(processed_dfs) * signs]
        centriods = centriods_array.mean(axis=0)
        centriod_df = pd.DataFrame(centriods,
                                   index=initial_df.index,
                                   columns=initial_df.columns)

        timer_print(calc_cent_start_time,
                    prefix="Centroid calculation is ready")
        # return the average components
        return centriod_df
