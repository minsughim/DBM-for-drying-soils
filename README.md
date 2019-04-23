# DBM-for-drying-soils
This repository includes the source codes for desert biocrust model under dynamic hydration conditions 


To excute the main function (Main_BSC_biogeoscience.m), you need four input arguments;

Time for the dynamics (examineDays)
matric potentail for unsaturated soils (pot1, with the unit of [-kPa])
time interval for output (plottt, with the unit of mins.)
index for ensemble averages (indexS, any integer to indicate a simulation) For example, to get the results of 5 days of simulation with output for every 5 mins at the matric potential of -3kPa, you excute in Matlab
examineDay = 5; 
pot1 = 3;
plottt = 5; 
indexS = 0;
Main_BSC_biogeoscience(examineDay, pot1, plottt, indexS)
Note

The main code include mex files that are complied for Mac and Linux systems.
