clear
load gridSiple1000
xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
[xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
uB = griddedInterpolant(xx,yy,tau(:,:,21));

tau_ISSM = uB(xy(:,1),xy(:,2));


old_tau = load('data/data_gridSipleXXSmall5000ISSM Shift2bedmap0.mat');
tau_c = old_tau.tau_c;
old_tau = tau_c(xy(:,1),xy(:,2),ones(size(xy(:,1))),ones(size(xy(:,1))));

spd = measures_interp('speed',xy(:,1),xy(:,2));
spd(isnan(spd)) = 1;

tau_shift = zeros(size(tau_ISSM));
tau_shift(spd<100) = tau_ISSM(spd<100).* 1.4; 
tau_shift(spd>100) = tau_ISSM(spd>100).* .7; 
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
tiledlayout(1,3)
nexttile(1)
% trisurf(t,xy(:,1),xy(:,2),uB(xy(:,1),xy(:,2)),'edgecolor','none')
trisurf(t,xy(:,1),xy(:,2),(tau_shift-old_tau)./uB(xy(:,1),xy(:,2)),'edgecolor','none')
colorbar
colormap(redblue)
title('Difference %')
view(2)
caxis([-1 1])
nexttile(2)
trisurf(t,xy(:,1),xy(:,2),tau_shift,'edgecolor','none')
colorbar
title('Shift')
view(2)
caxis([0 250e3])
nexttile(3)
trisurf(t,xy(:,1),xy(:,2),old_tau,'edgecolor','none')
colorbar
title('Old')
view(2)
caxis([0 250e3])



save('tauShiftable.mat','tau_ISSM','spd')
