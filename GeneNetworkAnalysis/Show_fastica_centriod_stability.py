"""
File:         Show_fastica_centriod_stability.py
Created:      2020-05-06
Last Changed: 2020-05-06
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
from classes.FastICA_wrapper import FastICA_wrapper
from classes.Centroid_decomposition_wrapper import Centroid_decomposition_wrapper
from classes.SVD_wrapper import SVD_types
import itertools
import numpy as np
import pandas as pd
import sys
import os
from utilities import read_pd_df, create_output_dir_if_not_exists
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    input_file_path = sys.argv[1]
    output_dir = sys.argv[2]
    create_output_dir_if_not_exists(output_dir)

    test_num_iterations = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
    n_runs_per_test = 25

    fastIca_Object = FastICA_wrapper(svd_type=SVD_types["FULL"],
                                     n_components=100)

    df = read_pd_df(input_file_path, {
                        "sep": "\t",
                        "index_col": 0
                    })
    comp_stabilities = {}
    for num_iterations in test_num_iterations:
        print("Perform runs with {} iterations".format(num_iterations))
        run_components = []
        for run_index in range(n_runs_per_test):
            print("Perform run {} of {}".format(run_index + 1, n_runs_per_test))
            centriod_wrapper = Centroid_decomposition_wrapper(fastIca_Object, n_iterations=num_iterations)
            centriod_wrapper.fit(df.T)

            output_components = centriod_wrapper.components
            run_components.append(output_components)
            comp_path = os.path.join(output_dir, "components_{}_{}.pickle".format(
                num_iterations, run_index + 1
            ))
            output_components.to_pickle(comp_path)
        comp_stabilities[num_iterations] = calc_comp_stability(run_components)

    stabilities = pd.DataFrame(comp_stabilities)
    stabilities.to_pickle(os.path.join(output_dir, "component_stabilities_df.pkl"))
    create_violinplot(stabilities, os.path.join(output_dir, "component_stabilities.jpg"))


def calc_comp_stability(runs, threshold=0.95):
    overall_ratios = []
    ica_funs_combinations = itertools.combinations(runs, 2)
    for run_comp_1, run_comp_2 in ica_funs_combinations:
        cor_mat = np.corrcoef(run_comp_1, run_comp_2)
        cor_mat = cor_mat[run_comp_1.shape[0]:, :run_comp_1.shape[0]]

        comp_above_threshold = np.abs(cor_mat).max(axis=0) >= threshold
        n_comp_above_treshold = np.sum(comp_above_threshold)
        ratio = n_comp_above_treshold / run_comp_1.shape[0]

        overall_ratios.append(ratio)
    return np.array(overall_ratios)

def create_violinplot(data, output_path=None):
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20, 10))
    sns.violinplot(data=data)
    plt.xlabel("Number of components")
    plt.ylabel("Ratio")
    plt.title(
        "Ratio of component pairwise correlations above threshold of 0.95 (25 fastica runs per number of components)")

    if output_path is not None:
        plt.savefig(output_path)
    else:
        plt.show()

if __name__ == '__main__':
    main()