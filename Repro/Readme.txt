These R-files and shell-script contain all code to reproduce the numerical results (figures/tables) for the paper 'Frequency Domain Statistical Inference for High-Dimensional Time Series'.

ALl necessary functions are listed in the R-package 'HDSpectralAnalysis', see HDSpectralAnalysis_0.1.0.9000.tar.gz or the github respiratory. Addtional packages are

abind 1.4-5
doRNG 1.8.6
dqrng 0.3.0
foreach 1.5.2
iterators 1.0.14
glmnet 4.1.7
glassoFast 1.0.1
mvtnorm 1.2.2
stats 4.3.1
network 1.14
R.matlab 3.7.0


Computation is done by running the shell-script (Simulations_ALL.sh). On a computer with 80 cores, this took about 2 weeks. Results are saved in RData files. 

The R-script Print_results_Paper_Repro.r reads the RData files and generates the figures and tables. 


For the real data example, first the data has to be downloaded at https://doi.org/10.18738/T8/CNVLAM and the subfolder Processed Empirical EEG Data: LAPref, Beta, Gaussian
The following files are used: 
EEG_Cat_Study4_Resting_Data_S19_EyesClosed_Gaussianized.mat
EEG_Cat_Study4_Resting_Data_S19_EyesOpen_Gaussianized.mat
EEG_Cat_Study4_Resting_Data_S20_EyesClosed_Gaussianized.mat
EEG_Cat_Study4_Resting_Data_S20_EyesOpen_Gaussianized.mat

The matlab files have to placed in a working-directory subfolder named Data_Subset. 

Then, the file Real_Data_Repro.r can be run. This file runs the analysis and produced the output in terms of latex-figures. 