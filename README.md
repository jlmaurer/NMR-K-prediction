# Predicting K from NMR data
This repo consists of code to find parameter values for models relating NMR measurements to hydraulic conductivity K. The models are given in Maurer and Knight, 2016. 

## Data
The data files are in the 'Aggregated_data' folder, and have columns: depth, T2ML, porosity, permeability, and optionally other measurements like Sum-of-Echoes (SOE). 

## Methods
Codes are provided for running grid search, bootstrap, MCMC, and direct solves for the parameters of the models. 

## Models
The codes are set up to handle the SDR model and various models using SOE. The Kozeny-Godefroy model can also be handled with some tweaking. 

## Licensing/Permissions/Acknowledgments
The data in the folder are the same as those published in Maurer and Knight 2016, see that paper for the original data sources and discussion. Codes were developed by and are copyrighed by Jeremy Maurer. All codes are freely available for use under the MIT license. 
