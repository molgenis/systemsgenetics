"""
File:         SVD_wrapper.py
Created:      2020-04-29
Last Changed: 2020-04-29
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
from sklearn.decomposition import PCA
from scipy import linalg
import numpy as np
import pandas as pd

SVD_types = {
    "AUTO": "auto",
    "FULL": "full",
    "RANDOM": "random",
    "SVD_GESDD": "svd_gesdd",
    "SVD_GESVD": "svd_gesvd"
}


class SVD_wrapper:
    def __init__(self, svd_type, n_components=None, white_data = False, output_dir = None):
        self.svd_type = svd_type
        self.n_components = n_components
        self.projected_data = None
        self.components = None
        self.explained_variance = None
        self.output_dir = output_dir
        self.white_data = white_data

    def perform_pca_auto(self, data):
        pca_object = PCA(n_components=self.n_components,
                         svd_solver="auto",
                         whiten=self.white_data)
        projected_data = pca_object.fit_transform(data)
        components = pca_object.components_
        self.explained_variance = pca_object.explained_variance_
        self.set_processed_results(data, components, projected_data)

    def perform_pca_full(self, data):
        pca_object = PCA(n_components=self.n_components,
                         svd_solver="full",
                         whiten=self.white_data)
        projected_data = pca_object.fit_transform(data)
        components = pca_object.components_
        self.explained_variance = pca_object.explained_variance_
        self.set_processed_results(data, components, projected_data)

    def perform_pca_random(self, data):
        pca_object = PCA(n_components=self.n_components,
                         svd_solver="randomized",
                         whiten = self.white_data,
                         iterated_power=20
                         )
        projected_data = pca_object.fit_transform(data)
        components = pca_object.components_
        self.explained_variance = pca_object.explained_variance_
        self.set_processed_results(data, components, projected_data)

    def perform_svd(self, data, lapack_driver="gesdd"):
        if self.n_components is None or self.n_components > data.shape[1]:
            self.n_components = data.shape[1]

        # Center the data
        centered_data = data - np.mean(data, axis=0)

        # perform SVD
        U, S, V = linalg.svd(centered_data,
                             full_matrices=False,
                             lapack_driver=lapack_driver)

        # flip eigenvectors' sign to enforce deterministic output
        max_abs_cols = np.argmax(np.abs(U), axis=0)
        signs = np.sign(U[max_abs_cols, range(U.shape[1])])
        U *= signs
        V *= signs[:, np.newaxis]

        # calculate the projected data
        projected_data = U[:, :self.n_components] * S[:self.n_components]
        components = V[:self.n_components]
        self.explained_variance = (S ** 2) / (data.shape[0] - 1)
        self.set_processed_results(data, components, projected_data)


    def set_processed_results(self, data, components, projected_data):
        components_index = pd.RangeIndex(start=1,
                                         stop=self.n_components + 1,
                                         name="IC")

        projected_data_df = pd.DataFrame(projected_data,
                                        index=data.index,
                                        columns=components_index)\
            .add_prefix("PC_")

        components_df = pd.DataFrame(components,
                                     index=components_index,
                                     columns=data.columns
                                     ).T.add_prefix("PC_").T

        self.components = components_df
        self.projected_data = projected_data_df


    def perform_svd_gesdd(self, data):
        self.perform_svd(data, lapack_driver="gesdd")

    def perform_svd_gesvd(self, data):
        self.perform_svd(data, lapack_driver="gesvd")

    def perform_run(self, data):
        if self.svd_type == SVD_types["AUTO"]:
            self.perform_pca_auto(data)
        if self.svd_type == SVD_types["FULL"]:
            self.perform_pca_full(data)
        if self.svd_type == SVD_types["RANDOM"]:
            self.perform_pca_random(data)
        if self.svd_type == SVD_types["SVD_GESDD"]:
            self.perform_svd_gesdd(data)
        if self.svd_type == SVD_types["SVD_GESVD"]:
            self.perform_svd_gesvd(data)

        if (self.output_dir is not None):
            np.savetxt(
                os.path.join(self.output_dir, 'explained_variance.csv'),
                self.explained_variance,
                delimiter=",")

    def fit(self, data):
        self.perform_run(data)

    def fit_transform(self, data):
        self.perform_run(data)
        return self.projected_data
