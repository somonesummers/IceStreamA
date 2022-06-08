clear; close all

xmax =  -0e5;
xmin = -8e5;
ymax =  -0e5;
ymin =  -8e5;

dx = 1e3;
smth = 5e3;
overgrab = 0;
xi = xmin-dx*overgrab:dx:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

b_raw =  bedmap2_interp(Xi,Yi,'bed');
sf_raw =  bedmap2_interp(Xi,Yi,'surface');
thickness =  bedmap2_interp(Xi,Yi,'thickness');
ice_mask =  bedmap2_interp(Xi,Yi,'icemask');

[u, v] = measures_interp('velocity',Xi,Yi);
spd = measures_interp('speed',Xi,Yi);
[lat, lon] = ps2ll(Xi,Yi);


% M = [Xi(:),Yi(:),lat(:),lon(:),u(:),v(:),spd(:),b_raw(:),sf_raw(:),thickness(:),ice_mask(:)];
% writematrix(M,'M.csv')


figure(6)
clf
ax1 = axes;
r = 8;
surf(ax1,Xi(1:r:end,1:r:end),Yi(1:r:end,1:r:end),b_raw(1:r:end,1:r:end),...
    ice_mask(1:r:end,1:r:end),'facealpha',[.9],'edgecolor', 'none');
% lighting gouraud
axis equal
c = colorbar('south');
% c.Label.String = 'Bed Elevation [m ASL]';
hold on
ax2 = axes;
p = surf(ax2,Xi,Yi,sf_raw,(spd),'facealpha',0.9);
title('Elevation with Bed')
set(p, 'edgecolor', 'none');
colormap(parula);
% caxis([-2000 2000])
axis equal
setFontSize(16);
c = colorbar('north');
c.Label.String = 'Ice Surface Speed [m/yr]';
f = gca;
f.ColorScale = 'log';
set(p, 'edgecolor', 'none');

hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'jet')
colormap(ax2,'parula')