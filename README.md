# DBM-for-drying-soils
![alt text](https://github.com/minsughim/DBM-for-drying-soils/blob/master/schematics_DBM_pH.PNG)

This repository includes source codes of the desert biocrust model (DBM) under dynamic hydration conditions.
This work is a continuation of the DBM at a static hydration condition [1]. You can find the main codes of this work at https://github.com/minsughim/DBM--Biogeoscience. To excute this model, follow the instruction there and excute the main code, Main_BSC_biogeoscience.m. 

The model for drying soils [2] uses the output of Main_BSC_biogeoscience.m as an initial condition. This ensures that distribtuion and activity of microbial communities are stabilised under given conditions (Think of sampling a real soil and performing experiments under desired conditions!). The stabilised communities undergo dry and dark conditions for 3 days (enough to exhaust carbon sources and nitrification activity becomes dominant). After the population dynamics reaches to a psuedo-steady state, a wet-dry cycle (duration of 24 hours) is applied. 

The main difference to the Main_BSC_biogeoscience.m is the dynamic update of hydration condition (assigned with matric potential). During wetting or drying, the domain conserves the masses of all compunds, leading to changes in concentrations and gases efflux. 

## Usage

### STEP 0. MATLAB Installation and download all files in the folder named 'codes'

### STEP 1. Inoculation of microbial cells at field capacity (-3kPa) 

Excute Main_BSC_biogeoscience.m and obtain the stabilised microbial community at field capacity.  
For this, you need four input arguments;
1. Time for the dynamics (examineDays)
2. matric potential for unsaturated soils (pot1, with the unit of [-kPa])
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
This code will generate a folder named "HT_BSC_Day5_Pot3.0_index0". 
Once the simulation is completed, there will be numtiple mat files in the folder, indluding:

Parameters.mat (physical domain, chemical parameters, and biological parameters will be saved here)
HTBioCrustCNcycleHour##.mat (intermediate results, that are saved with 12 hours interval:: only saves essential, population distribution, concentraiton of chemical compounds, etc.) 
Final.mat (The final result with all parameters to continue the simulation under different conditions:: about 1.5~2G depending on microbial population size)

### STEP 2. Microbial activity at fully saturated condtions

Once STEP1 is done, the simulation will continue for wetting of the domain by excuting following:
~~~~~~~~~~~~~{.m}
examineDay = 5; 
pot1 = 3; % at field capcity
plottt = 5; 
indexS = 0;

lightOn = 2; % 0 for complete darkness, 1 for constant light condition, any other numbers are for simulating diurnal cycles
examineDay2 = 5; %How many days to immerse the soil system in water?
Main_BSC_immerse(lightOn,examineDay,examineDay2 pot1, plottt, indexS)
~~~~~~~~~~~~~

### STEP 3, Incubation of the dry domain under darkness and apply a wetting-drying cycle
~~~~~~~~~~~~~{.m}
NH3ppb = 5; %atmospheric level of NH3 in ppb
HONOppb = 1; %atmospheric level of HONO in ppb
desiccationIndex = 4; %pre-calcuated drying patterns
newT = 25; %changing ambient temperature is also possible (degree celcius)
examineDay = 5; 
pot1 = 3; % at field capcity
indexS = 0;
Main_BSC_wet_dry(NH3ppb, HONOppb,desiccationIndex, newT, examineDay, pot1, indexS)
~~~~~~~~~~~~~

### STEP 4, Post processing of results 
~~~~~~~~~~~~~{.m}
hono = 1;
nh3 = 20;
desiccationIndex = 4;
newT = 25;
numberOfSamples = 8; % number of simulations under 1 boundary conditions
Process_results(hono, nh3,desiccationIndex,newT,numberOfSamples)
~~~~~~~~~~~~~


## Note

1. The DBM is computationally very expensive and rather slow on local desktops. For instance, using 32 cores in a computing cluster requires 4-5 days to complete a single simulation from STEP 1 to STEP 4

2. The main code includes mex files that are complied for Mac and Linux systems only (tested on MATLAB 2017b (9.3.0.713579) 64-bit)



## Reference

[1] Kim, M. and Or, D.: Hydration status and diurnal trophic interactions shape microbial community function in desert biocrusts, Biogeosciences, 14, 5403-5424, https://doi.org/10.5194/bg-14-5403-2017, 2017.

[2] Kim, M. and Or, D.: Microscale pH variations in drying soils and desert biocrusts affect HONO and NH3 emissions, in review, 2019.

