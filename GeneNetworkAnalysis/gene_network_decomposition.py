"""
File:         gene_network_decomposition.py
Created:      2020-04-30
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

import argparse
import logging
from utilities import create_output_dir_if_not_exists, create_logfile
from classes.Decomposition import Decomposition, decomposition_types
from classes.SVD_wrapper import SVD_types


def main():
    # parse the commandline arguments
    cmd_arguments = parse_arguments()

    # create output dir
    create_output_dir_if_not_exists(cmd_arguments.output)

    # activate logging
    create_logfile(cmd_arguments.output)

    try:
        decomp = Decomposition(
            input_file_path=cmd_arguments.input,
            output_dir=cmd_arguments.output,
            analysis_type=cmd_arguments.analysis,
            over_samples=cmd_arguments.samples,
            pca_type=cmd_arguments.pca_type,
            test_run=cmd_arguments.test,
            n_components=cmd_arguments.components,
            fastICA_max_iter=cmd_arguments.fastica_max_iter,
            fastICA_stable_iterations=cmd_arguments.fastica_stable_iterations,
            n_rows=cmd_arguments.nrows,
            perform_log2=cmd_arguments.log2,
            preprocessing_center_scale=cmd_arguments.center_scale,
            force=cmd_arguments.force,
            output_disable_gzip=cmd_arguments.output_disable_gzip,
            output_disable_txt=cmd_arguments.output_disable_txt,
            output_disable_pickle=cmd_arguments.output_disable_pickle
        )
        decomp.run()
    except Exception as ex:
        print(ex)
        logging.exception("Runtime error")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Perform the decomposition step of the updated '
                    'geneNetwork pipeline. The decomposition can be done by '
                    'using PCA (SVD with different algorithms), FastICA or '
                    'an stabilized version of FastICA')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        help='Input file path')

    parser.add_argument('-o', '--output',
                        required=True,
                        type=str,
                        help='Output directory path, '
                             'dir will be created if it doesn\'t exists')

    parser.add_argument('-a', '--analysis',
                        required=True,
                        type=str,
                        choices=decomposition_types.values(),
                        help='Select an analyse type')

    parser.add_argument('-p', '--q',
                        required=False,
                        type=str,
                        choices=SVD_types.values(),
                        help='Select the SVD type. The the analysis type is '
                             'set to FastICA, the selected PCA type will be '
                             'used to performing the whiting step.')

    parser.add_argument('-t', '--test',
                        required=False,
                        action='store_true',
                        help='Perform test run, use only 150 rows and 100 '
                             'features')

    parser.add_argument('-c', '--components',
                        required=False,
                        type=int,
                        help='Set the number of components')

    parser.add_argument('-m', '--fastica_max_iter',
                        required=False,
                        type=int,
                        help='Maximal number of FastICA iterations. '
                             'Can only be used if analysis is set to FastICA '
                             'or FastICA-stable. default is 2000')

    parser.add_argument('--fastica_stable_iterations',
                        required=False,
                        type=int,
                        help='Number of FastICA runs used to calculate the '
                             'final FastICA stable results. Can only be set '
                             'if analysis is set to FastICA-stable. '
                             'Default is 10')

    parser.add_argument('--fastica_stable_safe_intermediate',
                        required=False,
                        action='store_true',
                        help='Safe the intermediate FastICA runs in the '
                             'fastICA-stable run. Can only be set '
                             'if analysis is set to FastICA-stable.')

    parser.add_argument('-n', '--nrows',
                        required=False,
                        type=int,
                        help='Limit the number of rows in the input data')

    parser.add_argument('--input_transpose',
                        required=False,
                        action='store_true',
                        help='Transpose the input matrix before the analysis')

    parser.add_argument('-f', '--force',
                        required=False,
                        action='store_true',
                        help='Force to use the original input dataset instead '
                             'of cashed data')

    parser.add_argument('--log2',
                        required=False,
                        action='store_true',
                        help='Pre-processing input data with log2 '
                             'transformation')

    parser.add_argument('--center_scale',
                        required=False,
                        action='store_true',
                        help='Pre-processing input data with centering '
                             'and scaling')

    parser.add_argument('--output_disable_gzip',
                        required=False,
                        action='store_true',
                        help='Disable the gzip compression of the txt output files')

    parser.add_argument('--output_disable_txt',
                        required=False,
                        action='store_true',
                        help='Disable txt versions of the output files')

    parser.add_argument('--output_disable_pickle',
                        required=False,
                        action='store_true',
                        help='Disable pickled versions of the output files')

    return parser.parse_args()


if __name__ == '__main__':
    main()
