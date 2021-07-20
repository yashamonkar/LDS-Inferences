# LDS-Inferences

This is the code repository for the paper "A k-Nearest Neighbor Space-Time Simulator with applications to large-scale wind and solar power modelling" by Yash Amonkar (Columbia), David J Farnham (Carnegie Institue for Science) and Upmanu Lall (Columbia). 

All code to replicate the work done in the paper along with the supplementery materials is archived here. 
The scripts in the /codebase, when run, replicate the entire study along with additional cases and plots for an extensive analysis.
Most individual sub-routines are stored as functions in /functions.

The ERA-5 reanalysis data, for the Texas Interconnect region, is not added in this GitHub repository (due to size constraints) but is available at https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5. 
The original dataset is in  NetCDF (.nc) format and the NetCDF_to_Text.R script in data/rawdata is used to convert it to text (.txt) format. 
The entire datasets used and simulations generated is stored here (http://doi.org/10.5281/zenodo.5116072).
This is a version of the GitHub repository with data and simulations. Added separately due to the size constraints in GitHub


The four seperate R scripts in the codebase sub-folder perform the following functions.

1. KSTS_Simulator.R - Generates 48 independent simulations of the data using the novel KSTS simulator which also generates a number of plots analyzing the fit of the simulations to the data.

2. KNN_Simulator.R - Generates 48 independent simulations of the data using the novel KNN simulator. This also generates the same plots as above. 

3. Data_Characteristics.R - Explores the characteristcs of the reanalysis data. 

4. Simulation_Characteriscs.R - Explores the aggregated wind and solar profiles, Energy Droughts and Excedancee Probabilities. 

(Combined, theses four scripts generate all the analysis and plots used/displayed in the draft.)

This work has been submitted to Joule (https://www.cell.com/joule/home) for consideration.


Abstract:- 
The space-time variability of potential electricity generation  is a critical factor for the design of wind and solar dominated power systems. 
Limited historical data is available on incoming solar radiation and wind. 
Consequently, a stochastic model for simulating the co-variation of wind and solar fields is needed to assess the probability distribution of the severity and duration of  energy "droughts" at the network scale that need to be managed by long duration storage or alternate energy sources. 
We present a k-Nearest Neighbor Space-Time Simulator that accounts for the space-time dependence in high dimensional wind and solar fields and can simulate synthetic wind and solar fields of arbitrary length. 
Probability distributions of the severity-duration-frequency of regional energy droughts relative to a target are produced. 
As expected, the severity-duration-frequency of shortfalls is much greater when space-time dependence is properly considered than when it is not.
