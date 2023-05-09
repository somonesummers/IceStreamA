clear
load grids/gridFlowRiseC035.mat

overgrab = 20;
xxx = xmin-dx*overgrab:dx/2:xmax+dx*overgrab;
yyy = ymin-dx*overgrab:dx/2:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xxx,yyy);

rho   = 917;                % rho:   density of ice [kg/m^3]
rho_w = 1000;               % rho_w: density of water [kg/m^3]

% % % Raw fields
% bm2_b =  bedmap2_interp(Xi,Yi,'bed');
% bm2_s =  bedmap2_interp(Xi,Yi,'surface');
% 
% bm_b =  bedmachine_interp('bed',Xi,Yi);
% bm_s =  bedmachine_interp('surface',Xi,Yi);
% 
% [u,v] = measures_interp('velocity',Xi,Yi);
% spd2 = measures_interp('speed',Xi,Yi);



thin_m = 0;
str = 'ISSM';
mapFile = 'gridFlowRiseA02.mat';
initializeInputs();
initializeModel();
tau_c = defineTau(str);
buildSystem();

[u,v] = measures_interp('velocity',xy(:,1),xy(:,2));
u = u/3.154e7;
v = v/3.154e7;

figure(2)
clf
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'tight');

nexttile(1)

trisurf(t,xy(:,1),xy(:,2),h_bm_b(xy(:,1),xy(:,2)),...
               'edgecolor','none')
view(2)
colorbar
title('bed')

nexttile(2)
trisurf(t,xy(:,1),xy(:,2),h_bm_s(xy(:,1),xy(:,2)),...
               'edgecolor','none')
view(2)
colorbar
title('surf')

view(2)
colorbar
title('bed')

nexttile(3)
trisurf(t,xy(:,1),xy(:,2),h_bm(xy(:,1),xy(:,2)),...
               'edgecolor','none')
view(2)
colorbar
title('thickness')

nexttile(4)
trisurf(t,xy(:,1),xy(:,2),h_bm(xy(:,1),xy(:,2))+rho_w/rho*h_bm_b(xy(:,1),xy(:,2)),...
               'edgecolor','none')
view(2)
colorbar
title('overburden above floatation')



