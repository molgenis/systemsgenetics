"""
File:         FastICA_wrapper.py
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
from classes.SVD_wrapper import SVD_wrapper
from sklearn.decomposition import FastICA
from utilities import logging_print, stats_dict_to_string, timer_print
import numpy as np
import pandas as pd
import time


class FastICA_wrapper:

    def __init__(self, svd_type, n_components=None, max_iter=500):
        self.svd_type = svd_type
        self.n_components = n_components
        self.max_iter = max_iter

        self.whiten_components = None
        self.whiten_data = None

        self.projected_data = None
        self.components = None
        self.explained_variance = None

    def perform_data_whitening(self, data):
        pdw_start_time = time.time()
        logging_print("FastICA, start data whitening")
        whited_svd_object = SVD_wrapper(
            svd_type=self.svd_type,
            n_components=self.n_components,
            white_data=True
        )
        self.whiten_data = whited_svd_object.fit_transform(data).to_numpy()
        self.whiten_components = whited_svd_object.components.to_numpy()
        timer_print(pdw_start_time, prefix="FastICA, data whitening ready")

    def change_components_number(self, n_components):
        self.n_components = n_components
        if self.whiten_data.shape[1] <= n_components:
            self.whiten_components = None
            self.whiten_data = None

    def fit(self, data):
        if self.whiten_components is None:
            self.perform_data_whitening(data)
        fit_start_time = time.time()

        fastICA_object = FastICA(algorithm="parallel",
                                 whiten=False, fun='logcosh',
                                 max_iter=self.max_iter,
                                 tol=1e-10)

        fastICA_object.fit(self.whiten_data[:, :self.n_components])

        indep_comp = np.dot(fastICA_object.components_,
                            self.whiten_components[:self.n_components, :])
        indep_sources = np.dot(indep_comp, data.to_numpy().T).T

        components_index = pd.RangeIndex(start=1,
                                         stop=self.n_components + 1,
                                         name="IC")

        indep_sources_df = pd.DataFrame(indep_sources,
                                        index=data.index,
                                        columns=components_index
                                        ).add_prefix("IC_")

        indep_comp_df = pd.DataFrame(indep_comp,
                                     index=components_index,
                                     columns=data.columns
                                     ).T.add_prefix("IC_").T

        self.projected_data = indep_sources_df
        self.components = indep_comp_df

        logging_print(stats_dict_to_string({
            "Number of used iterations": fastICA_object.n_iter_
        }))
        timer_print(fit_start_time,
                    prefix="FastICA component optimalisation is ready")

    def fit_auto_white(self, data):
        logging_print("Use fastICA with auto whiten")
        fastICA_object = FastICA(n_components=self.n_components,
                                 algorithm="parallel",
                                 fun='logcosh',
                                 max_iter=500,
                                 tol=1e-10)

        sources = fastICA_object.fit_transform(data)
        self.projected_data = sources
        self.components = fastICA_object.components_

        logging_print(stats_dict_to_string({
            "Number of used iterations": fastICA_object.n_iter_
        }))

    def fit_transform(self, data):
        self.fit(data)
        return self.projected_data
