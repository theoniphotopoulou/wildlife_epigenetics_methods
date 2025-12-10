# wildlife_epigenetics_methods

Scripts used to produce results in "Methods for analysing wildlife DNA methylation data" by Photopoulou et al.

Main scripts: 

<!--- *00_data_prep.R* generates input data (mainly merges epigenetic data and health data and allocates fold IDs for CV) -->
- *00_make_synthetic_data.R* generates input data (creates and merges epigenetic data and health data and allocates fold IDs for CV) 
- *fgecs.R*, *probage.R*, to *healthmodels.R* can be run in standalone mode (need to remove comments from a few initial lines in code) or are called by *01_run_all.R*.
- *fgec.R* runs elastic net and random forest ("first generation") epigenetic clocks.
- *probage.R* runs a maximum likelihood implementation of Probage.
- *healthmodels.R* runs elastic net and random forest models with health outcomes rather than age as the response.
- *01_run_all.R* runs all modelling conditions. Runs over combinations of any hyperparameters so will take a long time - set up to run in parallel over 20-30 cores.
- *02_postprocessing.R* postprocesses results, primarily extracting test set results for optimized values of hyperparameters over the grid of values run in the previous step.
- *03_calc_predictions_and_aars.R* does additional postprocessing, calculating predicted ages and age acceleration residuals
- *04_plots_tables.R* produces tables and plots from the paper for epigenetic clock models.
- *05_plots_health.R* produces plots from the paper for trait score models.

Other scripts:

- *age-transformations-zoller.R*, *mammalmethylclockR_fns.R*, and *probage_fns.R* contain helper functions. *age-transformations-zoller.R* and *mammalmethylclockR_fns.R* are functions extracted from the **MammalMethylClock** R package, see [here](https://github.com/jazoller96/mammalian-methyl-clocks/tree/main/MammalMethylClock-Package).
- *run_existing_clocks.R* produces age estimates using functions from existing clocks using the functions stored in *mammalmethylclockR_fns.R*. These will be the same estimates obtained using the **MammalMethylClock** R package. Installing this package wasn't straightforward so we just extracted the necessary functions. Scripts require the *clocks_metadatabase.rda* file that contains clock coefficients, also extracted from **MammalMethylClock**.

