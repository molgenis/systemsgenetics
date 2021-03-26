"""
File:         create_auc_per_expl_var_plot.py
Created:      2020-12-17
Last Changed: 2020-12-17
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
import sys
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats

database_trans_table = {
"c2_cp_kegg_v7_1_entrez_gmt_matrix_txt": "KEGG",
"goa_human_2020_06_01_gaf_F_2020_06_01_matrix_txt": "GO",
"goa_human_2020_06_01_gaf_P_2020_06_01_matrix_txt": "GO",
"phenotype_to_genes_V1268_OMIMandORPHA_txt_matrix": "HPO",
"Ensembl2Reactome_All_Levels_2020_07_18_txt_matrix_txt": "Reactome",
"goa_human_2020_06_01_gaf_C_2020_06_01_matrix_txt": "GO",
}

color_table = {
    "KEGG": "#3AA976",
    "GO": "#1275B0",
    "HPO": "#F7941D",
    "Reactome": "#AF14A4",
}

last_row_index = None

def calculate_p_val(total_df):
    groups = total_df.loc[:, "group"].unique()
    print("groups", groups)
    last_index = None
    for row_index in groups:
        if last_index != None:
            group_1_values = total_df.loc[total_df["group"] == last_index, "AUC"]
            group_2_values = total_df.loc[total_df["group"] == row_index, "AUC"]

            t_stat, p_val = scipy.stats.ttest_ind(group_1_values.dropna(), group_2_values.dropna())
            print("{} against {}\t{}".format(last_index, row_index, p_val))
        last_index = row_index


def main():
    expl_var_dirs = glob.glob(os.path.join(sys.argv[1], "*"))
    auc_tables = []

    for expl_var_path in expl_var_dirs:
        if os.path.isdir(expl_var_path):
            expl_variance = expl_var_path.split("_")[4][1]
            expl_variance_title = "0.{}".format(expl_variance)
            print(expl_var_path, expl_variance_title)
            auc_file_paths = glob.glob(os.path.join(expl_var_path, "*", "predictions_auc_bonf.pkl"))
            for auc_file_path in auc_file_paths:
                auc_df = pd.read_pickle(auc_file_path)
                database_name = auc_file_path.split("/")[2]
                database_title = database_trans_table[database_name]
                print(database_title, database_name)
                single_df = pd.DataFrame(
                    {"AUC": auc_df["auc"],
                     "pathway": database_title,
                     "group": expl_variance_title})
                auc_tables.append(single_df)


    total_df = pd.concat(auc_tables)
    total_df.sort_values(by="group", inplace=True)
    fig, ax = plt.subplots(1, 1, figsize=(1 * 14, 12), subplot_kw= dict(facecolor="white"))

    calculate_p_val(total_df)

    print("total_df means average", total_df.groupby(["group"]).mean())
    sns.lineplot("group", "AUC", hue="pathway", style="pathway",
                 markers=["o"] * 4, dashes=False,
                 data=total_df,
                 palette=color_table,
                 # ci="sd",
                 ci=None,
                 sort=False,
                 zorder=10, ax=ax, alpha=0.3, err_kws={"alpha": 0.3})
    sns.lineplot("group", "AUC", markers=["o"] * 4, dashes=False,
                 data=total_df, sort=False, ci=None,
                 zorder=10, ax=ax, color="black", label="Average")
    ax.set_xlabel("Explained variance")
    ax.set_ylabel("Mean AUC")
    ax.set_title("Explained variance cutoffs of GeneNetwork data")

    ax.grid(False)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))
    ax.spines['left'].set_color("black")
    ax.spines['bottom'].set_color("black")
    plt.legend(title="Pathway name")




    fig_path = os.path.join(sys.argv[1], "genenetwork_expl_var_auc.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print("Figure created: {}".format(fig_path))
    fig_path = os.path.join(sys.argv[1], "genenetwork_expl_var_auc.pdf")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print("Figure created: {}".format(fig_path))

if __name__ == '__main__':
    main()