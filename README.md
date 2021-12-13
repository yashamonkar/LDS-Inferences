# LDS-Inferences

This is the code repository for the paper "A k-Nearest Neighbor Space-Time Simulator with applications to large-scale wind and solar power modelling" by Yash Amonkar (Columbia), David J Farnham (Carnegie Institue for Science/ClimateAi) and Upmanu Lall (Columbia). 
This work has been submitted to Patterns (https://www.cell.com/patterns/home) for consideration and is in the advanced revision/pre-publication stage.


All code to replicate the work done in the paper along with the supplementery materials is archived here. 
The scripts in the /codebase, when run, replicate the entire study along with additional cases and plots for an extensive analysis.
Most individual sub-routines are stored as functions in /functions.

The ERA-5 reanalysis data, for the Texas Interconnect region, is not added in this GitHub repository (due to size constraints) but is available at https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5. 
The original dataset is in  NetCDF (.nc) format and the NetCDF_to_Text.R script in data/rawdata is used to convert it to text (.txt) format. 
The entire datasets used and simulations generated are stored here - (http://doi.org/10.5281/zenodo.5116072).
This is a version of the GitHub repository with data and simulations. Added separately due to the size constraints in GitHub


The four seperate R scripts in the codebase sub-folder perform the following functions.

1. KSTS_Simulator.R - Generates 48 independent simulations of the data using the novel KSTS simulator which also generates a number of plots analyzing the fit of the simulations to the data.

2. KNN_Simulator.R - Generates 48 independent simulations of the data using the novel KNN simulator. This also generates the same plots as above. 

3. Data_Characteristics.R - Explores the characteristcs of the reanalysis data. 

4. Simulation_Characteriscs.R - Explores the aggregated wind and solar profiles, Energy Droughts and Excedancee Probabilities. 

(Combined, theses four scripts generate all the analysis and plots used/displayed in the draft.)


Abstract:- 
We develop and present a k-Nearest Neighbor Space-Time Simulator that accounts for the spatio-temporal dependence in high dimensional hydroclimatic fields (e.g. wind and solar) and can simulate synthetic realizations of arbitrary length.
We illustrate how this statistical simulation tool can be used in the context of regional power system planning under a scenario of high reliance on wind and solar generation and when long historical records of wind and solar power generation potential are not available.
We show how our simulation model can be used to assess the probability distribution of the severity and duration of energy ``droughts" at the network scale that need to be managed by long duration storage or alternate energy sources.
We present this estimation of supply side shortages for the Texas Interconnection.