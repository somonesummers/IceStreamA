% Plot figures of region for overview figures
clear; close all
xbox = [-4.9500  -3.6000   -2.7000   -4.1000   -4.9500]*1e5;
ybox = [-3.8000  -3.5300   -5.2250   -5.6700   -3.8000]*1e5;

xmax =  -2.55e5;
xmin = -5.05e5;
ymax =  -3.75e5;
ymin =  -5.75e5;

pv = [xbox; ybox]';

icey = cbrewer('div','BrBG',48);
rho = 917;
rho_w = 1000;
g = 9.81;
B = 1.6e8; % A = 2.4e-25 Pa^(-3) s^(-1)
A = 3.5e-25;
overgrab = 0;
% xmax =  -2.5e5;
% xmin =  -5.0e5;
% ymax =  -3.5e5;
% ymin =  -6.0e5;

dx = 500;
% smth = 6e3;
xi = xmin-dx*overgrab:dx:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% Masking
% msk = scatteredInterpolant(xy(:,1),xy(:,2),ones(size(xy(:,1))),'nearest','none');
% mask = msk(Xi,Yi);
% mask(isnan(mask)) = 0;
 mask = ones(size(Xi)); %make mask all ones for now

%% Raw fields
[u, v] = measures_interp('velocity',Xi,Yi);
b_raw =  bedmachine_interp('bed',Xi,Yi);
sf_raw =  bedmachine_interp('surface',Xi,Yi);

%% Cleaning NANs in velocity [do cautiously]
valid_u     = ~isnan(u);
u_interp = scatteredInterpolant(Xi(valid_u),Yi(valid_u),u(valid_u),'natural');
u(~valid_u) = u_interp(Xi(~valid_u),Yi(~valid_u));
valid_v     = ~isnan(v);
v_interp = scatteredInterpolant(Xi(valid_v),Yi(valid_v),v(valid_v),'natural');
v(~valid_v) = v_interp(Xi(~valid_v),Yi(~valid_v));

%% Smoothing
sf_smooth = imgaussfilt(sf_raw,5e3/dx);
b_smooth = imgaussfilt(b_raw,5e3/dx);
u_smooth = imgaussfilt(u,2e3/dx);
v_smooth = imgaussfilt(v,2e3/dx);
% b_raw =  bedmap2_interp(Xi,Yi,'bed');
% sf_raw =  bedmap2_interp(Xi,Yi,'surface');

h = sf_smooth - b_smooth;
[sx ,  sy] = gradient(sf_smooth,dx,dx);
[ux ,  uy] = gradient(u_smooth,dx,dx);
[vx ,  vy] = gradient(v_smooth,dx,dx);
e_eff = sqrt(.5*(ux.^2 + vy.^2) + (.5*(ux + uy).^2) + (.25*(uy + vx)).^2);
spd = sqrt(u.^2 + v.^2);
alpha = atan(v./u);
e_xy = .5 * (uy + vx);
e_shr = (vy-ux).*cos(alpha).*sin(alpha) + e_xy.*(cos(alpha).^2 - (sin(alpha).^2));
Tdx = -rho * g * h .* sx;
Tdy = -rho * g * h .* sy;
Td  = sqrt(Tdx.^2 +  Tdy.^2);
e_xz = A * Td.^3 * 3.154e7; %[1/year]

sf_ridge = zeros(size(sf_raw));
sf_ridge(spd < 30) = sf_raw(spd < 30);
%%

figure('Position',[50 500 800 600])

surf(Xi/1e3,Yi/1e3,zeros(size(spd)),spd,'edgecolor','none');
hold on
axis equal
c = colorbar;
c.Label.String = 'Ice Surface Speed [m/yr]';
caxis([1 700])
ax = gca;
ax.ColorScale = 'log';
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
[cc, hh] = contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
% % clabel(cc,hh,'LabelSpacing',500,'Color','k','Fontsize',20)
% contour(xi/1e3,yi/1e3,sf_raw, (0:50:700) , '-','HandleVisibility','off','linewidth',2,'color',rgb('gray'))
% contour(xi/1e3,yi/1e3,sf_raw, (0:10:700) , '-','HandleVisibility','off','color',rgb('gray'))
view(2)
xlabel('Easting [km]')
ylabel('Northing [km]')
setFontSize(28)
mapzoomps('nw','km')
% savePng('figs/RiseAMap')
contour(xi/1e3,yi/1e3,sf_raw, (0:50:700) , '-','HandleVisibility','off','linewidth',2,'color',rgb('gray'))
contour(xi/1e3,yi/1e3,sf_raw, (0:10:700) , '-','HandleVisibility','off','color',rgb('gray'))
% savePng('figs/RiseAMapContour')

%%
figure('Position',[50 500 800 600])

surf(Xi/1e3,Yi/1e3,40*ones(size(e_shr)),e_shr,'edgecolor','none','facealpha',.8);
hold on
colormap(cbrewer('div','RdBu',128))
c = colorbar;
axis equal
c.Label.String = 'Lateral Shear Strain Rate [1/yr]';
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
[cc, hh] = contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
view(2)
caxis([-0.03 0.03])
setFontSize(28)
xlabel('Easting [km]')
ylabel('Northing [km]')

mapzoomps('nw','km')
% savePng('figs/RiseAMapShear')

%% 
figure('Position',[50 500 800 600])
H_ab = (sf_raw - b_raw) + 1000/917*b_raw;
surf(Xi/1e3,Yi/1e3,H_ab,'edgecolor','none','facealpha',.8);
hold on
colormap(cbrewer('seq','OrRd',128));
c = colorbar;
axis equal
axis off
c.Label.String = 'Height Above Floatation [m]';
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
[cc, hh] = contour(xi/1e3,yi/1e3,H_ab, [50:50:250], '-','linewidth',2,'HandleVisibility','off','color',rgb('light red'));
clabel(cc,hh,'LabelSpacing',1000,'Color','k','Fontsize',20)

contour(xi/1e3,yi/1e3,H_ab, [0, 0], '-','linewidth',2,'HandleVisibility','off','color',rgb('black'));
view(2)
caxis([0 500])
setFontSize(28)
xlabel('Easting [km]')
ylabel('Northing [km]')
savePng('figs/RiseAMapHAF')
% mapzoomps('nw','km')

%%
figure('Position',[50 500 800 600])

surf(Xi/1e3,Yi/1e3,zeros(size(spd)),spd,'edgecolor','none');
hold on
axis equal
c = colorbar;
c.Label.String = 'Ice Surface Speed [m/yr]';
caxis([1 700])
ax = gca;
ax.ColorScale = 'log';
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
% contour(xi/1e3,yi/1e3,spd, [30:10:700] , '-','linewidth',1,'HandleVisibility','off','color',rgb('light gray'));
% contour(xi/1e3,yi/1e3,spd, [100:100:700] , '-','linewidth',2,'HandleVisibility','off','color',rgb('light gray'));
grid1 = load("grids/gridFlowRiseA02.mat");
plot(grid1.pv(:,1)/1e3,grid1.pv(:,2)/1e3,'r--','linewidth',2)
view(2)
xlabel('Easting [km]')
ylabel('Northing [km]')
setFontSize(28)
% mapzoomps('sw','km')
% savePng('figs/domainMap')

%%
% tau_c = defineTau('ISSM');
% 
% figure('Position',[50 500 800 600])
% 
% surf(Xi/1e3,Yi/1e3,zeros(size(spd)),tau_c(Xi,Yi,1,1)/1e3,'edgecolor','none');
% shading flat
% hold on
% axis equal
% colormap(cmocean('thermal'))
% c = colorbar;
% c.Label.String = 'Basal Strength [kPa]';
% caxis([1 150])
% ax = gca;
% % ax.ColorScale = 'log';
% contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
% contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
% % contour(xi/1e3,yi/1e3,spd, [30:10:700] , '-','linewidth',1,'HandleVisibility','off','color',rgb('light gray'));
% % contour(xi/1e3,yi/1e3,spd, [100:100:700] , '-','linewidth',2,'HandleVisibility','off','color',rgb('light gray'));
% % grid1 = load("grids/gridFlowRiseA02.mat");
% % plot(grid1.pv(:,1)/1e3,grid1.pv(:,2)/1e3,'r--','linewidth',2)
% view(2)
% xlabel('Easting [km]')
% ylabel('Northing [km]')
% setFontSize(28)
% % savePng('figs/taucMap')

%%
figure('Position',[50 500 1600 600])
tiledlayout(1,3, 'Padding', 'tight', 'TileSpacing', 'tight');
clf
sgtitle('Strain Estimates')
ax = nexttile(1);
surf(Xi/1e3,Yi/1e3,40*ones(size(e_xz)),e_xz,'edgecolor','none','facealpha',.9);
hold on
colormap(ax,cbrewer('seq','Reds',128))
% colormap(ax,'parula')
c = colorbar;
axis equal
title('Vertical Shear Strain Rate')
c.Label.String = '[1/yr]';
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
[cc, hh] = contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
view(2)
clim([2e-5 2e-3]);
setFontSize(14)
xlabel('Easting [km]')
ylabel('Northing [km]')

ax = nexttile(2);
surf(Xi/1e3,Yi/1e3,40*ones(size(e_eff)),abs(e_eff),'edgecolor','none','facealpha',.9);
hold on
colormap(ax,cbrewer('seq','Reds',128))
% colormap(ax,'parula');
c = colorbar;
axis equal
c.Label.String = '[1/yr]';
title('Horizontal Effective Strain Rate')
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
[cc, hh] = contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
view(2)
clim([2e-5 2e-3]);
setFontSize(14)
xlabel('Easting [km]')
ylabel('Northing [km]')



ax = nexttile(3);
surf(Xi/1e3,Yi/1e3,40*ones(size(e_xz)),e_eff./e_xz,'edgecolor','none','facealpha',.8);
hold on
tempMap = (cbrewer('seq','Oranges',400));
tempMap2 = flipud(cbrewer('seq','Purples',50));
colormap(ax,[tempMap2;tempMap])
% colormap(ax,'parula')
c = colorbar;
axis equal
title("Ratio of Effective Horizontal Strain to Vertical Strain");
c.Label.String = "[ ]";
contour(xi/1e3,yi/1e3,spd, (0:2:30) , 'k-','HandleVisibility','off')
[cc, hh] = contour(xi/1e3,yi/1e3,spd, [10,20,30] , 'k-','linewidth',2,'HandleVisibility','off');
view(2)
clim([0 10])
setFontSize(14)
xlabel('Easting [km]')
ylabel('Northing [km]')

% mapzoomps('nw','km')
% savePng('figs/RiseAMapShear')