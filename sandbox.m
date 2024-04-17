clear
close all

file = "data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0";
data1 = load("data/spdChange/" + file + ".mat");

figure
ax = gca;
plotLogStrain(data1,3,ax)