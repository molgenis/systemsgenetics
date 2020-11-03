General structure:

I generally use .rmd files for downstream interpretation of results, while .r files are scripts to be run on the cluster.
Generally I like being able to use the code blocks r markdown offers.

In the data folder there are source_x.r files. Here you can put the paths to the relevant files on your specific machine.
This is to avoid constantly overwriting eachothers paths in the git.

I'm not putting the .RPoject and .Rhistory files in the git. Personally I have these in the downstreamer_main folder.


└── downstreamer_main
    ├── analysis                  > Folder for a given analysis
    │   ├── data            > Small data files and scripts to source data for specific user
    │   └── output          > Folder for output, output files not per se in git
    │       └── plots       > Folder for plots, output files not per se in git


Current tree:

├── downstreamer_collaborations     > Scripts for collaborators
│   ├── biogen
│   ├── covid_consortium
│   ├── john_perry
│   └── lude_cancer_gwas
└── downstreamer_main               > Scripts for main analysis
    ├── evaluating_coregulation     > Coregulation related plots / analyses
    ├── gwas_simulations            > Scripts for doing the GWAS simulations, and evaluating the results
    ├── legacy_scripts              > Legacy scripts not yet organized or deprecated
    │   ├── examples
    │   │   ├── coeliac
    │   │   ├── height
    │   │   └── kidney
    │   └── testing_scripts
    ├── network_plots               > Scripts for making networks & plots
    └── umap                        > Scripts relating to making uMAP