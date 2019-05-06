clear all
maxT = 44;
poolobj = parpool('local',maxT); 
HONOppb = 0.005;
NH3ppb = 20;
Nitrate_singlePatch2(NH3ppb, HONOppb)