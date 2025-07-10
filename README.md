# wildlife_epigenetics_methods

Scripts to reproduce results in "Methods for analysing wildlife DNA methylation data" by Photopoulou et al.

Main scripts: 

- *00_data_prep.R* generates input data (mainly merges epigenetic data and health data and allocates fold IDs for CV)
- *fgecs.R*, *probage.R*, to *healthmodels.R* can be run in standalone mode (need to remove comments from a few initial lines in code) or are called by *01_run_all.R*.
- *fgec.R* runs elastic net and random forest ("first generation") epigenetic clocks.
- *probage.R* runs a maximum likelihood implementation of Probage.
- *healthmodels.R* runs elastic net and random forest models with health outcomes rather than age as the response.
- *01_run_all.R* runs all modelling conditions. Runs over combinations of any hyperparameters so will take a long time - set up to run in parallel over 20-30 cores.
- *02_postprocessing.R* postprocesses results, primarily extracting test set results for optimized values of hyperparameters over the grid of values run in the previous step.
- *03_plots_tables.R* produces tables and plots from the paper (and some others not used).

Other scripts:

- *age-transformations-zoller.R* and *probage_fns.R* contain helper functions.
- *run_existing_clocks.R* produces age estimates for existing clocks using the **MammalMethylClock** R package. Installing this package wasn't straightforward and needed various bioconductor packages so suggest just using results in *output/existing_clocks.Rdata* (these are used only in *03_plots_tables.R*).

