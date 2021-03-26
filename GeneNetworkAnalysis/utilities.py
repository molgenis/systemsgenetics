"""
File:         utilities.py
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
import time
from datetime import datetime
import logging
import pandas as pd
import numpy as np


def read_pd_df(file_path, pd_read_csv_args=None, force=False, no_error=False,
               proc_df_before_save=None):
    if os.path.isfile(file_path):
        pickle_file_path = "{}_cashed.pickle".format(file_path)
        if os.path.isfile(pickle_file_path) and force is False:
            return pd.read_pickle(pickle_file_path)
        else:
            if file_path.endswith(".pkl") or file_path.endswith(".pickle"):
                data = pd.read_pickle(file_path)
                if proc_df_before_save is not None:
                    data = proc_df_before_save(data)
                return data
            else:
                if pd_read_csv_args is None:
                    pd_read_csv_args = {}
                data = pd.read_csv(file_path, **pd_read_csv_args)
                if proc_df_before_save is not None:
                    data = proc_df_before_save(data)
                data.to_pickle(pickle_file_path)
                return data
    else:
        if no_error == False:
            raise FileNotFoundError("Cannot find input file: {}".format(
                file_path
            ))
    return None


def stats_dict_to_string(stats_info):
    input_string = []
    for key, value in stats_info.items():
        new_value = value
        if isinstance(value, list):
            new_value = ", ".join(value)
        input_string.append("{}: {}".format(key, new_value))

    return "\n".join(input_string)


def logging_print(text_line, print_line=True, log_line=True):
    if print_line:
        print(text_line)
    if log_line:
        logging.info(text_line)


def timer_print(start_time,
                end_time=None, prefix=None, print_info=True,
                time_overview_log=None):
    if end_time is None:
        end_time = time.time()
    if prefix is None:
        prefix = "Ready"
    proc_min, proc_sec = divmod(end_time - start_time, 60)
    line = "{prefix} (in {min:0.0f} min and " \
           "{sec:0.0f} sec)".format(prefix=prefix,
                                    min=proc_min,
                                    sec=proc_sec)
    logging_print(line, print_line=print_info)
    if time_overview_log is not None and prefix is not None:
        proc_time_in_sec = proc_min * 60  + proc_sec
        if prefix in time_overview_log:
            if isinstance(time_overview_log[prefix], list):
                time_overview_log[prefix].append(proc_time_in_sec)
            else:
                time_overview_log[prefix] = [time_overview_log[prefix],
                                             proc_time_in_sec]
        else:
            time_overview_log[prefix] = proc_time_in_sec


def create_output_dir_if_not_exists(dir_path):
    if os.path.isdir(dir_path) is False:
        os.mkdir(dir_path)
    return dir_path


def create_logfile(output_dir):
    log_file_path = os.path.join(output_dir, "analysis_log.log")
    logging.basicConfig(filename=log_file_path,
                        filemode='a',
                        level=logging.DEBUG,
                        format='%(asctime)s %(message)s')


def log_matrix_info(df, name, log_info=True, print_info=True):
    df_info = "{name} n_row: {n_row}, n_col: {n_col}\n" \
              "first column headers: {column_index}\n" \
              "first row index: {row_index}" \
              "".format(name=name,
                        n_row=df.shape[0],
                        n_col=df.shape[1],
                        column_index=df.columns.values[:5],
                        row_index=df.index.values[:5]
                        )
    logging_print(df_info, print_line=print_info, log_line=log_info)


def remove_file_if_present(file_path):
    try:
        if file_path is not None and os.path.isfile(file_path):
            os.remove(file_path)
    except OSError:
        logging_print("cannot remove file: '{}'".format(file_path))


def remove_gene_from_pathway(hpo_term_list, gene):
  hpo_list_adj = hpo_term_list.copy()
  hpo_list_adj[gene] = False
  return hpo_list_adj


def create_file_path(dir_name, file_name, extension=None):
    ex_file_name = file_name
    ex_file_name = ex_file_name.format(datetime=datetime.now())
    if extension is not None:
        ex_file_name = "{}.{}".format(ex_file_name, extension)
    return os.path.join(dir_name, ex_file_name)


def create_auc_output_file(auc_values, gene_pathway_count, p_values,
                                output_dir, bonf_p_values=None,
                                file_name=None, export_txt=True, txt_gzip=True,
                                export_pkl=True):
    predictions_auc_table = pd.DataFrame({
        "geneCount": gene_pathway_count,
        "pValue": p_values,
        "auc": auc_values,
    })
    predictions_auc_table.index.name = "term"

    if bonf_p_values is not None:
        predictions_auc_table["bonferroni"] = bonf_p_values

    export_file_name = file_name
    if file_name is None:
        if bonf_p_values is not None:
            export_file_name = "predictions_auc_bonf"
        else:
            export_file_name = "predictions_auc"

    if export_txt:
        if txt_gzip:
            predictions_auc_table.to_csv(
                create_file_path(output_dir, export_file_name, 'txt.gz'),
                sep='\t')
        else:
            predictions_auc_table.to_csv(
                create_file_path(output_dir, export_file_name, 'txt'),
                sep='\t')
    if export_pkl:
        predictions_auc_table.to_pickle(
            create_file_path(output_dir, export_file_name, 'pkl'))


def create_gene_pathway_output_file(gene_pathway_scores, output_dir,
                                    file_name=None, export_txt=True,
                                    txt_gzip=True, export_pkl=True):
    export_file_name = file_name
    if file_name is None:
        export_file_name = "gene_pathway_scores"

    if export_txt:
        if txt_gzip:
            gene_pathway_scores.to_csv(
                create_file_path(output_dir, export_file_name, "txt.gz"),
                sep='\t')
        else:
            gene_pathway_scores.to_csv(
                create_file_path(output_dir, export_file_name, "txt"),
                sep='\t')

    if export_pkl:
        gene_pathway_scores.to_pickle(
            create_file_path(output_dir, export_file_name, "pkl"))


def remove_gene_from_hpo(hpo_term_list, gene):
  hpo_list_adj = hpo_term_list.copy()
  hpo_list_adj[gene] = False
  return hpo_list_adj


def create_output_dir_name(output_dir_name, multi_node_node_id=None,
                           multi_node_num_nodes=None):
    if multi_node_node_id is not None and multi_node_num_nodes is not None:
        create_output_dir_if_not_exists(output_dir_name)
        return os.path.join(output_dir_name,
                            "node_{}".format(multi_node_node_id+1))
    return output_dir_name