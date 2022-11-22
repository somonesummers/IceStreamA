clear
load grids/gridRefinedRiseF05.mat

overgrab = 20;
xi = xmin-dx*overgrab:dx/2:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx/2:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% % % Raw fields
% bm2_b =  bedmap2_interp(Xi,Yi,'bed');
% bm2_s =  bedmap2_interp(Xi,Yi,'surface');
% 
% bm_b =  bedmachine_interp('bed',Xi,Yi);
% bm_s =  bedmachine_interp('surface',Xi,Yi);

[u,v] = measures_interp('velocity',Xi,Yi);
spd = measures_interp('speed',Xi,Yi);

% pv = makePerimeter(pv, 10e3);

x   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
y   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
[xx,yy] = ndgrid(x - 3072000,y - 3072000);
tic
uA = griddedInterpolant(xx,yy,tau(:,:,21));
a = uA(0,0);
toc
tic
uB = griddedInterpolant({x - 3072000,y - 3072000},tau(:,:,21));
b = uB({0,0});
toc


%%
figure(1)
clf
quiver(Xi,Yi,u,v)
hold on
contour(xi,yi,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title('Quiver')
plot(pv(:,1),pv(:,2),'-*b')
h1 = streamline(xi,yi,u,v,pv(1,1)-1e4,pv(1,2));
h2 = streamline(xi,yi,u,v,pv(2,1),pv(2,2));
set(h1,'Color','none');
set(h2,'Color','none');
plot(h1.XData(1:100:2700),h1.YData(1:100:2700),'*-','color',rgb('Forest Green'))
plot(h2.XData(1:100:1300),h2.YData(1:100:1300),'*-','color',rgb('Forest Green'))
pv_new = [h1.XData(1:100:2700)',h1.YData(1:100:2700)';flipud(h2.XData(1:100:1300)'),flipud(h2.YData(1:100:1300)');h1.XData(1), h1.YData(1)];
plot(pv_new(:,1),pv_new(:,2))
% h = streamline(xi,yi,-u,-v,pv(:,1),pv(:,2));
% set(h,'Color',rgb('purple'));
% h = streamline(xi,yi,-u,-v,-2.95e5:.5e3:-2.9e5,-5.15e5:.5e3:-5.1e5);
% set(h,'Color',rgb('lavender'));
axis equal
save('pv_new','pv_new')
