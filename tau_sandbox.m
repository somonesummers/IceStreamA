clear
load gridSiple1000
xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
[xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
uB = griddedInterpolant(xx,yy,tau(:,:,21));

tau_ISSM = uB(xy(:,1),xy(:,2));

spd = measures_interp('speed',xy(:,1),xy(:,2));
spd(isnan(spd)) = 1;

tau_shift = zeros(size(tau_ISSM));
tau_shift(spd<100) = tau_ISSM(spd<100).* 1.7; 
tau_shift(spd>100) = tau_ISSM(spd>100).* .3; 
figure(10)
clf
scatter(spd,tau_ISSM,'filled')
hold on
scatter(spd,tau_shift,'filled')
% plot(spd,,'k--')
xlabel('speed')
ylabel('strength')

figure(11)
clf
tiledlayout(1,2)
nexttile(1)
trisurf(t,xy(:,1),xy(:,2),uB(xy(:,1),xy(:,2)),'edgecolor','none')
colorbar
title('ISSM Inversion')
view(2)
caxis([0 250e3])
nexttile(2)
trisurf(t,xy(:,1),xy(:,2),tau_shift,'edgecolor','none')
colorbar
title('Shift')
view(2)
caxis([0 250e3])



save('tau_shift2.mat','tau_shift')
