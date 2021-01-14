"""
File:         gene_pathway_score_prediction.py
Created:      2020-02-28
Last Changed: 2021-01-08
Author(s):    H.H. Wiersma

Copyright (C) 2019 H.H. Wiersma

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.`

A copy of the GNU General Public License can be found in the LICENSE file in
the root directory of this source tree. If not, see
<https://www.gnu.org/licenses/>.
"""


import os
import logging
import argparse
from classes.Gene_pathway_analysis import Gene_Pathway_Analysis, analysis_types
from utilities import create_output_dir_if_not_exists, create_output_dir_name


def main():
    cmd_arguments = parse_arguments()

    if "," in cmd_arguments.pathway_matrix:
        pathway_paths = cmd_arguments.pathway_matrix.split(",")
        if cmd_arguments.genes is not None and "," in cmd_arguments.genes:
            pathway_genes_paths = cmd_arguments.genes.split(",")
            if len(pathway_paths) == len(pathway_genes_paths):
                if cmd_arguments.node is not None and cmd_arguments.num_nodes is not None:
                    nodes_per_pathway = cmd_arguments.num_nodes // len(pathway_paths)
                    pathway_id, node_id = divmod(cmd_arguments.node, nodes_per_pathway)
                    pathway_path = pathway_paths[pathway_id]
                    pathway_genes = pathway_genes_paths[pathway_id]
                    pathway_name = "_".join(os.path.basename(pathway_path).split(".")[:-1])
                    create_output_dir_if_not_exists(cmd_arguments.output)
                    base_dir = os.path.join(cmd_arguments.output, pathway_name)

                    # create output dir
                    output_dir = create_output_dir_name(
                        base_dir,
                        node_id,
                        nodes_per_pathway)
                    create_output_dir_if_not_exists(output_dir)

                    # activate logging
                    log_file_path = os.path.join(output_dir,
                                                 "analysis_log.log")
                    logging.basicConfig(filename=log_file_path,
                                        filemode='a',
                                        level=logging.DEBUG,
                                        format='%(asctime)s %(message)s')
                    try:
                        analysis = Gene_Pathway_Analysis(
                            components_data_path=cmd_arguments.input,
                            matrix_path=pathway_path,
                            output_dir=output_dir,
                            analysis_type=cmd_arguments.analysis,
                            minimal_number_of_genes=cmd_arguments.minimal_genes,
                            background_genes_path=pathway_genes,
                            n_cores=cmd_arguments.cores,
                            split_start=cmd_arguments.start,
                            split_end=cmd_arguments.end,
                            multi_node_node_id=node_id,
                            multi_node_num_nodes=nodes_per_pathway,
                            multi_node_output_dir=base_dir,
                            force=cmd_arguments.force,
                            output_disable_gzip=cmd_arguments.output_disable_gzip,
                            output_disable_txt=cmd_arguments.output_disable_txt,
                            output_disable_pickle=cmd_arguments.output_disable_pickle
                        )
                        analysis.run()
                    except Exception as ex:
                        print(ex)
                        logging.exception("Runtime error")
    else:
        # create output dir
        output_dir = create_output_dir_name(cmd_arguments.output,
                                            cmd_arguments.node,
                                            cmd_arguments.num_nodes)
        create_output_dir_if_not_exists(output_dir)

        # activate logging
        log_file_path = os.path.join(output_dir, "analysis_log.log")
        logging.basicConfig(filename=log_file_path,
                            filemode='a',
                            level=logging.DEBUG,
                            format='%(asctime)s %(message)s')
        try:
            analysis = Gene_Pathway_Analysis(
                components_data_path=cmd_arguments.input,
                matrix_path=cmd_arguments.pathway_matrix,
                output_dir=output_dir,
                analysis_type=cmd_arguments.analysis,
                minimal_number_of_genes=cmd_arguments.minimal_genes,
                background_genes_path=cmd_arguments.genes,
                n_cores=cmd_arguments.cores,
                split_start=cmd_arguments.start,
                split_end=cmd_arguments.end,
                multi_node_node_id=cmd_arguments.node,
                multi_node_num_nodes=cmd_arguments.num_nodes,
                multi_node_output_dir=cmd_arguments.output,
                force=cmd_arguments.force,
                output_disable_gzip=cmd_arguments.output_disable_gzip,
                output_disable_txt=cmd_arguments.output_disable_txt,
                output_disable_pickle=cmd_arguments.output_disable_pickle
            )
            analysis.run()
        except Exception as ex:
            print(ex)
            logging.exception("Runtime error")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Perform a fastica or pca over the gene sample matrix '
                    'in the gene network pipeline')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        help='Components input file')

    parser.add_argument('-o', '--output',
                        required=True,
                        type=str,
                        help='Path to the output directory')

    parser.add_argument('-p', '--pathway_matrix',
                        required=True,
                        type=str,
                        help='Path to the pathway matrix file')

    parser.add_argument('-g', '--genes',
                        required=False,
                        type=str,
                        help='File with the background genes')

    parser.add_argument('-a', '--analysis',
                        required=True,
                        type=str,
                        choices=analysis_types.values(),
                        help='Select a analyse type')

    parser.add_argument('-m', '--minimal_genes',
                        required=False,
                        type=int,
                        help='Minimal number of genes in the matrix annotation')

    parser.add_argument('--start',
                        required=False,
                        type=int,
                        help='First id of the pathways to start analysing')

    parser.add_argument('--end',
                        required=False,
                        type=int,
                        help='Last id of the pathways to analyse')

    parser.add_argument('--node',
                        required=False,
                        type=int,
                        help='Node id to split the analysis over multiple nodes')

    parser.add_argument('--num_nodes',
                        required=False,
                        type=int,
                        help='Total number of nodes')

    parser.add_argument('-c', '--cores',
                        required=False,
                        type=int,
                        help='Number of cores for multiprocessing')

    parser.add_argument('-f', '--force',
                        required=False,
                        action='store_true',
                        help='Force to use the orginal input dataset instead of pickled data')

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
