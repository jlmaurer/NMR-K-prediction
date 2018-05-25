%% Script for testing different sites using SOE
%% Load Data
clear

name = 'dpnmr_leque_all'; 

% load data file
[d, Dk, ~, phi, z, SumEch, kk, ~, lp, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(name); 
logSumEch = log10(SumEch); 

