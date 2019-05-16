# DBM-for-drying-soils
This repository includes source codes of the desert biocrust model (DBM) under dynamic hydration conditions.
This work is a continuation of the DBM at a static hydration condition [1]. You can find the main codes of this work at https://github.com/minsughim/DBM--Biogeoscience. To excute this model, follow the instruction there and excute the main code, Main_BSC_biogeoscience.m. 

The model for drying soils [2] uses the output of Main_BSC_biogeoscience.m as an initial condition. This ensures that distribtuion and activity of microbial communities are stabilised under given conditions (Think of sampling a real soil and performing experiments under desired conditions!). This stabilised communities undergo dry and dark conditions for 3 days (enough to exhaust carbon sources and nitrification activity becomes dominant). After the population dynamics reaches to a psuedo-steady state, a wet-dry cycle (duration of 24 hours) is applied. 

The main difference to the Main_BSC_biogeoscience.m is the dynamic update of the hydration conditions of the domain. During wetting or drying, the domain conserves the masses of all compunds, leading to changes in concentrations and gases efflux. 

## Usage

STEP 1. Inoculation of microbial cells at field capacity (-3kPa) 
Excute Main_BSC_biogeoscience.m and obtain the stabilised microbial community.  
For this, you need four input arguments;
1. Time for the dynamics (examineDays)
2. matric potentail for unsaturated soils (pot1, with the unit of [-kPa])
3. time interval for output (plottt, with the unit of mins.)
4. index for ensemble averages (indexS, any integer to indicate a simulation)
For example, to get the results of 5 days of simulation with output for every 5 mins at the matric potential of -3kPa, you excute in Matlab
~~~~~~~~~~~~~{.m}
examineDay = 5; 
pot1 = 3; % at field capcity
plottt = 5; 
indexS = 0;
Main_BSC_biogeoscience(examineDay, pot1, plottt, indexS)
~~~~~~~~~~~~~
STEP 2. Microbial activity at fully saturated condtions


STEP 3, 'Incubation' of the system under dark and dry condition


STEP 4, Apply a wetting-drying cycle



## Note

1. The DBM is computationally very expensive and rather slow on local desktops. For instance, using 32 cores in a computing cluster requires 4-5 days to complete a single simulation from STEP 1 to STEP4 with a following procedure

2. The main code include mex files that are complied for Mac and Linux systems (tested on MATLAB 2018a)

## Reference

[1] Kim, M. and Or, D.: Hydration status and diurnal trophic interactions shape microbial community function in desert biocrusts, Biogeosciences, 14, 5403-5424, https://doi.org/10.5194/bg-14-5403-2017, 2017.

[2] Kim, M. and Or, D.: Microscale pH variations in soils affect HONO and NH3 emissions during drying, in review, 2019.

