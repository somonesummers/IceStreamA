clear
close all
load gridSiple5000.mat

xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
[xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
a = load('bedDragDx5000smth10000Capped.mat');

aB = griddedInterpolant(a.Xi',a.Yi',subplus(a.bed'));

figure
trisurf(t,xy(:,1),xy(:,2),uB(xy(:,1),xy(:,2)),'edgecolor','none')
colorbar
title('ISSM Inversion')
view(2)
caxis([0 250e3])


figure
trisurf(t,xy(:,1),xy(:,2),aB(xy(:,1),xy(:,2)),'edgecolor','none')
colorbar
title('VDV Method')
view(2)
caxis([0 250e3])

figure
trisurf(t,xy(:,1),xy(:,2),aB(xy(:,1),xy(:,2))-uB(xy(:,1),xy(:,2)),'edgecolor','none')
colorbar
title('VDV - ISSM')
colormap redblue
view(2)
caxis([-100e3 100e3])