clear
close all
load('/Users/paulsummers/Documents/IceStreamA/grids/gridFlowRiseA02.mat')

figure(1)
% trisurf(t,xy(:,1)/1e3,xy(:,2)/1e3,ep(xy(:,1),xy(:,2)))
spd_new = fillmissing(measures_interp('speed',xy(:,1),xy(:,2)),'linear');
trisurf(t,xy(:,1)/1e3,xy(:,2)/1e3,spd_new,'facealpha',.7,'facecolor','interp')
view(2)
xlabel('Easting [km]')
ylabel('Northing [km]')
setFontSize(18)
% mapzoomps('nw','km')
% colormap(cbrewer('seq','YlGn',128))
c = colorbar;
c.Label.String = 'Observed Surface Speed [m/yr]';
axis equal
savePng('figs/domainAndSpeed')

