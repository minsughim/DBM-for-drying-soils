# DBM-for-drying-soils
This repository includes the source codes for desert biocrust model under dynamic hydration conditions.
The work is a continuation work of the DBM at a static hydration condition [1], (the main code for this work is in other repository can be run with Main_BSC_biogeoscience.m). 
The model uses the output of Main_BSC_biogeoscience.m as an initial condition, used as the stabilised microbial communities under the given conditions, and the domain is under the dry and dark condition (for 3 days; enough to exhaust carbon sources and nitrification activity becomes dominant), then a wet-dry cycle (duration of 24 hours) is applied. 




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

1. The DBM is computationally very expensive and rather slow on local desktops. For instance, using 32 cores in a computing cluster requires 4-5 days to complete a single simulation with a following procedure; (1) inoculation of microbial cells at field capacity (-3kPa), (2) calculation of the pseudo-steady state of microbial activities at fully saturated conditions, (3) simulation of a dry condition under darkness, and (4) application of a wetting-drying cycle.

2. The main code include mex files that are complied for Mac and Linux systems (tested on MATLAB 2018a)


## Reference

[1] Kim, M. and Or, D.: Hydration status and diurnal trophic interactions shape microbial community function in desert biocrusts, Biogeosciences, 14, 5403-5424, https://doi.org/10.5194/bg-14-5403-2017, 2017.

[2] Kim, M. and Or, D.: Microscale pH variations in soils affect HONO and NH3 emissions during drying, in review, 2019.

