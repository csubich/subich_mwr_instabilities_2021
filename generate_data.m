% Script to re-generate data for all growth-rate figures.  Takes about 24h
% to run per figure.

clear
fourier_instab
save growthrate_fourier_512
clear
linear_instab
save growthrate_linear_512
clear
cubic_instab
save growthrate_subic_512