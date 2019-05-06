pause(30)
clear all
maxT = 44;
poolobj = parpool('local',maxT);
HONOppb = 1;
NH3ppb = 5;
Nitrate_singlePatch2(NH3ppb, HONOppb)
