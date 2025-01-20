clear
close all
load('/Users/paulsummers/Documents/IceStreamA/grids/gridFlowRiseA02.mat')

xmin = min(xy(:,1));
xmax = max(xy(:,1));
ymin = min(xy(:,2));
ymax = max(xy(:,2));

dx = 1e3;
over = 25;
x = (xmin - over*dx):dx:(xmax + over*dx);
y = (ymin - over*dx):dx:(ymax + over*dx);
[xx,yy] = meshgrid(x,y);
spd = measures_interp('speed',xx,yy);


figure('Position',[908 395 1071 772])
surf(x/1e3,y/1e3,zeros(size(spd)),spd,'edgecolor','none')
hold on
% trisurf(t,xy(:,1)/1e3,xy(:,2)/1e3,ep(xy(:,1),xy(:,2)))
% spd_new = fillmissing(measures_interp('speed',xy(:,1),xy(:,2)),'linear');
trisurf(t,xy(:,1)/1e3,xy(:,2)/1e3,ones(size(xy(:,1))),'edgealpha',.7,'facecolor','none')
% plot3(pv(:,1)/1e3,pv(:,2)/1e3,2*ones(size(pv(:,1))),'r--','linewidth',3)
contour3(x/1e3,y/1e3,spd,[5,15,25],'k')
contour3(x/1e3,y/1e3,spd,[10,20,30],'k','linewidth',2)


view(2)
xlabel('Easting [km]')
ylabel('Northing [km]')
setFontSize(18)
% mapzoomps('nw','km')
% colormap(cbrewer('seq','YlGn',128))
c = colorbar;
c.Label.String = 'Observed Surface Speed [m/yr]';
ax = gca;
ax.ColorScale = 'log';
xlim([xmin-over*dx xmax+over*dx]/1e3);
ylim([ymin-over*dx ymax+over*dx]/1e3);
axis equal
title('Model Domain and Grid')
savePng('figs/domainAndSpeed')

