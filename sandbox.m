clear
close all

file = "data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0";
load("data/spdChange/" + file + ".mat");

fg1 = figure;
fg2 = figure;

inLoopPlotting

lambda = calcAdvection(T,u,v,xy,dx,rho,C_p);
La =@(x,y) lambda(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y)));

figure
trisurf(t,xy(:,1),xy(:,2),La(xy(:,1),xy(:,2)),'edgecolor','none')
view(2)
colorbar
title('surf T')

T_adj = max(T,T_bar(xy(:,1),xy(:,2))); % Ice can't be colder than surf
lambda = calcAdvection(T_bar(xy(:,1),xy(:,2)),u,v,xy,dx,rho,C_p);
La =@(x,y) lambda(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y)));
 
figure
trisurf(t,xy(:,1),xy(:,2),La(xy(:,1),xy(:,2)),'edgecolor','none')
view(2)
colorbar
title('T bar')

figure
trisurf(t,xy(:,1),xy(:,2),T_bar(xy(:,1),xy(:,2)),'edgecolor','none')
view(2)
colorbar

figure
trisurf(t,xy(:,1),xy(:,2),T,'edgecolor','none')
view(2)
colorbar