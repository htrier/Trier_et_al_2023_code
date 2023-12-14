This directory contains code and data to reproduce results presented in the manuscript titled ‘An ancient subcortical circuit decides when to orient to threat in humans’ (Trier, Khalighinejad, Hamilton, Harbison, Priestley, Laubach, Scholl, & Rushworth, 2023).

## System requirements
-	This code was developed on Apple silicon M2 Mac with Ventura 13.5 OS.
- R version 4.3.2 arm64
- RStudio version 2023.09.1+494
- R packages: `rjson 4.3.0`, `data.table 4.3.11`, `dplyr 4.3.1`, `tibble 3.2.1`, `tictoc 1.2`, `brms 2.20.4`, `ggplot2 3.4.4`, `rstatix 0.7.2`, `ggpubr 0.6.0`, `plotrix 3.8.4`, `emmeans 1.8.9`, `reshape 0.8.9`.
- Python version 3.11.6
- Python packages: `IPython 8.16.1`, `ipykernel 6.25.2`, `jupyter_client 8.3.1`, `jupyter_core 5.3.2`, `jupyter_server 2.7.3`, `jupyterlab 4.0.6`, `nbclient 0.8.0`, `nbconvert 7.9.2`, `nbformat 5.9.2`, `notebook 7.0.4`, `traitlets 5.10.1`, `seaborn : 0.13.0`, `pandas : 2.1.1`, `numpy: 1. 26.0`, `scipy: 1.11.3`, `matplotlib : 3.8.2`, `scikit_posthocs : 0.8.0`, `scikit-learn : 1.3.1`, `statsmodels : 0.14.0`
- Matlab R2023b with Statistics and Machine Learning Toolbox

## Installation guide
- Download R here: https://cran.r-project.org) and RStudio here: https://posit.co/download/rstudio-desktop/. 
- Install the required R packages by entering install.packages(‘name_of_package’) in the RStudio console.
- Install Jupyter by entering ‘pip3 install jupyter’ at the command line.
- Install each python package by entering ‘pip3 install package==version’ at the command line.
- Download Matlab here: https://uk.mathworks.com/downloads/web_downloads/?s_tid=hp_ff_t_downloads

## Instructions to visualize and reproduce results
Instructions are organized by the figure in which the result is shown. Full fMRI data cannot be included in this format, but we have included in the `/data/` folder materials demonstrating the output of the GLM model applied to the fMRI data:
- `ROI_stats_all_participants.csv`, which contains ROI statistics from the whole-brain analysis extracted at the individual level. This is the same information as Table S6 in the supplementary materials.
- `all-ROI-stats-2Mar2023-indexed-thresh0.0001_formatted.xlsx`, which contains ROI statistics from the whole-brain analysis at the group level.
- `all-cluster-stats-thresh0.0001_2Feb2023.csv`, which contains all cluster statistics at the p<0.0001 threshold from the whole-brain analysis at the group level.
- zstat files at the group level for contrasts of interest are included in `/Trier_et_al_2023_code/data/whole_brain_GLM_group_stats/`, in folders titled by contrast.

### Fig. 1H-L
-	Open `processRawBehavioralData.Rmd` in RStudio. Update the root paths for your computer. Run each cell to read the raw behavioural data from the experimental session (i.e. button presses and task contextual information), compute variables (e.g., behavioural variables such as which check is the first in a sequence, and computational variables such as timePressure), and export the preprocessed data to csv files. Not all behavioural measures will be used in the analysis. Output from this script is already saved in `/data/`. Runtime: 240s.
-	Open a terminal window, navigate to the `/Trier_et_al_2023_code/` folder, and type ‘Jupyter notebook’. Then open the notebook analyse_RTs.ipynb. Update the root paths for your computer. This notebook reads in the preprocessed behavioral data from each participant creates one large dataframe with all reaction times organized by the type of action switch they represent (forage->forage, forage->check, check->check, check->forage). It also creates dataframes checkSeqRTs and forageSeqRTs which log the reaction times with respect to the position of each button press within a sequence of actions. These data frames are saved for you in the /data/ folder and will be used in the next analysis step. Runtime: 1h.
-	Open `model_behavior.Rmd`. Update the root paths for your computer. Run all cells to model choice and RT data and reproduce figures 1H-L. Runtime: When re-running all bayesian models, runtime is up to 10 hours depending on system specs. However, the model output is shared in the `/data/` folder and the code has an option to read in that data and quickly reproduce the figures with it.
-	Reproducibility: Please note that Markov chain Monte Carlo is a stochastic algorithm. In theory each time the algorithm is run it generates a different Markov chain and therefore different Markov chain Monte Carlo estimators. These estimators will converge to the same answer only for infinitely long Markov chains. Thus, finite Markov chains and the subsequent estimators are not fully reproducible. We have included in the `/results/` folder the model results reported in our manuscript. Re-running the code should produce a very similar result. In addition, the R function ggplot does not plot jittered data points at the same exact position every time, so plots may look slightly different than the manuscript. However, testing our model results with the code will reproduce the ttest statistics reported in the supplementary materials.

### Figs. 2D-E, 6C-D
-	Open a terminal window, navigate to the `/Trier_et_al_2023_code/` folder, and type `jupyter notebook`. Then open the notebook `Plot_and_test_avg_ROI_activations.ipynb`. Update the root paths for your computer. Running all cells in the notebook will reproduce the bar charts from Fig. 2D-E, 6C-D, and create a statistics table called ‘all_stats’.

### Figs. 3a-b, 7Aiii-Aiv
-	To reproduce the correlations shown in Figures 3a, 3b, 7Aiii, and 7Aiv, open a terminal window, nativate to the `/Trier_et_al_2023_code/` folder, and type `jupyter notebook`. Then open the notebook `Figs_3a_3b_7Aiii_7Aiv.ipynb`. Update the root paths for your computer. Then execute each cell and the figures and corresponding statistics tables will be produced as output. Some text in the notebook further describes the behavioural measures being used.

### ANOVAs shown in Figs. 4Ai, 4Bi, 5A, 5E
-	To reproduce the ANOVAs shown in Figs. 4Ai, 4Bi, 5A, 5E, open a terminal window, navigate to the `/Trier_et_al_2023_code/` folder, and type `jupyter notebook`. Then open the notebook `ANOVAs_Figs_4Ai_4Bi_5A_5E.ipynb`. Update the root path for your computer. Then execute each cell and tables showing the corresponding ANOVAs and post-hoc tests will be produced as output.

### Figs. 4Aii-iv, 4Bii-iv, 7Ai-Aii, 5b-d, 5f-k, 7Bi-Bii, 7C, 7D
To reproduce the bar charts from panels 7Ai-Aii, 7Bi-Bii, and 7C, open each of the following scripts in Matlab. Update the root paths for your computer. Ensure that the folder `/Trier_et_al_2023_code/` is included in your Matlab path because it contains the data folder and supporting Matlab functions. Ensure that all .zip files in the `/data/` folder are unzipped. Then each script should be able to execute and produce the desired figure without further changes.
- Fig_4Aii.m
- Fig_4Aiii.m
- Fig_4Aiv.m
- Fig_4Bii.m
- Fig_4Biii.m
- Fig_4Biv.m
- Fig_5bcd.m
- Fig_5fgh.m
- Fig_5ijk.m
- Figs_7Ai_7Aii.m
- Fig_7Bi.m
- Fig_7Bii.m
- Fig_7C.m
- Fig_7D.m
