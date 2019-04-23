# DBM-for-drying-soils
This repository includes the source codes for desert biocrust model under dynamic hydration conditions 

## Usage

To excute the main function (Main_BSC_biogeoscience.m), you need four input arguments;
1. Time for the dynamics (examineDays)
2. matric potentail for unsaturated soils (pot1, with the unit of [-kPa])
3. time interval for output (plottt, with the unit of mins.)
4. index for ensemble averages (indexS, any integer to indicate a simulation)
For example, to get the results of 5 days of simulation with output for every 5 mins at the matric potential of -3kPa, you excute in Matlab

~~~~~~~~~~~~~{.m}
examineDay = 5; 
pot1 = 3;
plottt = 5; 
indexS = 0;
Main_BSC_biogeoscience(examineDay, pot1, plottt, indexS)
~~~~~~~~~~~~~

## Note

The main code include mex files that are complied for Mac and Linux systems.


## Reference

[1] Kim, M. and Or, D.: Hydration status and diurnal trophic interactions shape microbial community function in desert biocrusts, Biogeosciences, 14, 5403-5424, https://doi.org/10.5194/bg-14-5403-2017, 2017.
[1] Kim, M. and Or, D.: Microscale pH variations in soils affect HONO and NH3 emissions while drying, in review, 2019.

