"""
File:         cutoff_eigenvectors.py
Created:      2020-07-08
Last Changed: 2020-07-08
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

import argparse
import numpy as np
import pandas as pd


def main():
    cmd_arguments = parse_arguments()
    eig_cutter = EigenvectorsCutter(
        input_file_path=cmd_arguments.input,
        output_file_path=cmd_arguments.output,
        components_number=cmd_arguments.number,
        variance_cutoff=cmd_arguments.cutoff,
        explained_variance_file_path=cmd_arguments.variance
    )
    eig_cutter.run()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Script to limit the eigenvectors file to a given number '
                    'of components (--number) or to a given explained '
                    'variance ratio (--cutoff)')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        help='Path to eigenvectors file (pickled file or txt')

    parser.add_argument('-v', '--variance',
                        required=False,
                        type=str,
                        help='Path to explained variance file')

    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        help='Path to the output eigenvectors file with the '
                             'limited number of components. If no path is '
                             'given, the orginal file name will be used '
                             'with suffix "_filtered"')
    parser.add_argument('-n', '--number',
                        required=False,
                        type=int,
                        help='Cutoff of number of components')
    parser.add_argument('-c', '--cutoff',
                        required=False,
                        type=float,
                        help='Cutoff value ratio of explained variance of the '
                             'components. value between 0 and 1')
    return parser.parse_args()


class EigenvectorsCutter:
    def __init__(self, input_file_path, output_file_path=None,
                 components_number=None, variance_cutoff=None,
                 explained_variance_file_path=None):

        self.input_file_path = input_file_path
        self.output_file_path = output_file_path
        self.components_number = components_number
        self.variance_cutoff = variance_cutoff
        self.explained_variance_file_path = explained_variance_file_path

        self.input_df = None

    def run(self):
        self.read_input_file()
        self.filter_dataframe()
        self.create_output_files()

    def read_input_file(self):
        extension = self.input_file_path.split(".")[-1]
        if extension in ["pkl", "pickle"]:
            self.input_df = pd.read_pickle(self.input_file_path)
        else:
            self.input_df = pd.read_csv(self.input_file_path, sep='\t')

    def get_export_file_path(self, extension="pkl"):
        if self.output_file_path is not None and self.output_file_path != "":
            return "{}.{}".format(self.output_file_path, extension)

        export_path_name = ".".join(self.input_file_path.split(".")[:-1])
        return "{}__filtered.{}".format(export_path_name, extension)

    def filter_dataframe(self):
        if self.variance_cutoff is not None and self.variance_cutoff != "" and \
                self.explained_variance_file_path is not None and \
                self.explained_variance_file_path != "":
            explain_var_df = pd.read_csv(self.explained_variance_file_path,
                                         sep=",",
                                         header=None,
                                         names=["explVar"])

            total_variance = np.sum(explain_var_df.iloc[:, 0], axis=None)
            cunsum_expl_variance = np.cumsum(
                explain_var_df.iloc[:, 0]) / total_variance

            indx_of_component = np.argmax(cunsum_expl_variance >=
                                          self.variance_cutoff)

            print("filtered based on {:0.02f} total explaind variance. "
                  "Component number {} will be "
                  "used as cutoff".format(self.variance_cutoff,
                                          indx_of_component))

            self.input_df = self.input_df.iloc[:, :indx_of_component]
        elif self.components_number is not None and \
                self.components_number != "":
            print("filtered on component number: {}".format(
                self.components_number))
            self.input_df = self.input_df.iloc[:, :self.components_number]

    def create_output_files(self):
        self.input_df.to_csv(self.get_export_file_path("txt"),
                             sep='\t')
        self.input_df.to_pickle(self.get_export_file_path("pkl"))
        print("export files created: {}  and {}".format(
            self.get_export_file_path("txt"),
            self.get_export_file_path("pkl")
        ))


if __name__ == '__main__':
    main()
