clear
close all

file = "data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0";
data1 = load("data/spdChange/" + file + ".mat");

[u,v] = measures_interp('velocity',data1.xy(:,1),data1.xy(:,2));

data1.u = u/3.154e7;
data1.v = v/3.154e7;

xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];

figure
tiledlayout(1,3)
ax = nexttile(1);
plotLogStrain(data1,1,ax)
xlim(xLimits)
ylim(yLimits)
title('Long')

ax = nexttile(2);
plotLogStrain(data1,2,ax)
xlim(xLimits)
ylim(yLimits)
title('Trans')

ax = nexttile(3);
plotLogStrain(data1,3,ax)
xlim(xLimits)
ylim(yLimits)
title('Shear')

setFontSize(24)



