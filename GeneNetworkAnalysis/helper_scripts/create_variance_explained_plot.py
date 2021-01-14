"""
File:         create_variance_explained_plot.py
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
import matplotlib.pyplot as plt


def main():
    cmd_arguments = parse_arguments()

    variance_plot_maker = VariancePlot(
        input_path=cmd_arguments.input,
        output_path=cmd_arguments.output,
        max_components=cmd_arguments.max,
        visualize_cutoff=cmd_arguments.cutoff
    )
    variance_plot_maker.run()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Create a variance explained plot from a given csv file '
                    'from the decomposition script')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        help='Path to input file')
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        help='Path to the output file. If no path is given, the'
                             ' image will get the same name as the input file')
    parser.add_argument('-m', '--max',
                        required=False,
                        type=int,
                        help='Maximal number of components to plot')
    parser.add_argument('-c', '--cutoff',
                        required=False,
                        type=int,
                        help='Number of a particular component to highlight in'
                             ' the plot (for illustrating a cutoff value)')
    return parser.parse_args()


class VariancePlot:
    def __init__(self, input_path, output_path=None, max_components=None,
                 visualize_cutoff=None):
        self.input_path = input_path
        self.output_path = output_path
        self.max_components = max_components
        self.visualize_cutoff = visualize_cutoff
        self.input_df = None

    def get_export_img_path(self):
        if self.output_path is not None and self.output_path != "":
            return self.output_path
        export_path = ".".join(self.input_path.split(".")[:-1])
        return "{}.png".format(export_path)

    def run(self):
        self.read_input()
        self.filter_on_max_number_of_components()
        self.create_plot()

    def read_input(self):
        self.input_df = pd.read_csv(self.input_path,
                                    sep=",",
                                    header=None,
                                    names=["explVar"])

        self.total_variance = np.sum(self.input_df.iloc[:, 0], axis=None)

    def filter_on_max_number_of_components(self):
        if self.max_components is not None and self.max_components != "":
            self.input_df = self.input_df.iloc[:self.max_components, :]

    def create_plot(self):
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))

        cunsum_expl_variance = np.cumsum(self.input_df.iloc[:, 0]) / \
                               self.total_variance

        ax.plot(range(1, cunsum_expl_variance.shape[0] +1),
                cunsum_expl_variance, '*-', c="lightseagreen", zorder=1)

        if self.visualize_cutoff is not None and \
                self.visualize_cutoff != "" and \
                self.visualize_cutoff < self.input_df.shape[0]:
            var_expl_by_cutoff = cunsum_expl_variance[self.visualize_cutoff]

            ax.scatter(self.visualize_cutoff +1, var_expl_by_cutoff,
                       s=100, c='tomato', zorder=10)
            ax.axhline(var_expl_by_cutoff,
                       c="black", dashes=(5, 5), alpha=0.5, zorder=5)

            print("Total explained variance by cutoff "
                  "(component: {}) is: {:0.02f}".format(self.visualize_cutoff,
                                                        var_expl_by_cutoff))
        ax.set_xlabel("Components")
        ax.set_ylabel("Explained variance ratio")
        ax.set_title("Cumulative explained variance")
        export_img_path = self.get_export_img_path()
        plt.savefig(export_img_path, dpi=300, bbox_inches='tight')
        print("Explained variance plot is created, path: {}".format(
            export_img_path))


if __name__ == '__main__':
    main()
