# LDS-Inferences

This is the code repository for the paper "Inferences on Long-Duration Storage Requirements for Renewable Energy Systems from a Novel Space-Time Simulator" by Yash Amonkar (Columbia), David J Farnham (Carnegie Institue for Science) and Upmanu Lall (Columbia).

The ERA-5 reanalysis data, for the Texas Interconnect region, is not added in this GitHub repository but is available at https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5. The original dataset is in  NetCDF (.nc) format and the NetCDF_to_Text.R script in data/rawdata is used to convert it to text (.txt) format. 

This work has been submitted to Joule (https://www.cell.com/joule/home) for consideration. 


Abstract:- 
The space-time variability of potential electricity generation is a critical factor for the
design of wind and solar dominated power systems. Limited historical data is available
on incoming solar radiation and wind. Consequently, a stochastic model for simulating
the co-variation of wind and solar fields is needed to assess the probability distribution
of the severity and duration of energy "droughts" at the network scale that need to be
managed by long duration storage or alternate energy sources. We present a novel k-
Nearest Neighbor Space-Time Simulator that accounts for the space-time dependence
in high dimensional wind and solar fields and can simulate synthetic wind and solar
fields of arbitrary length. Probability distributions of the severity-duration-frequency of
regional energy droughts relative to a target are produced. As expected, the severityduration-
frequency of shortfalls is much greater when space-time dependence is
properly considered than when it is not.
