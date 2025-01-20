%% Plot all figures for 1D change, This is best plotting option
clear
close all
addpath lib/
saveFigs = true;

% if(saveFigs)
%     disp("Please confirm you'd like to save figures");
%     pause()
% end

[Acc, T_s] = loadALBMAP();


%% Cases of thickness
groupName = 'DhDt_';
% cases = [5:-1:-1]; %speed cases
baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";

cValues = [-3 3];
cMapToUse = flipud(cbrewer('div','PiYG',256));
data2 = load(baseFile);

%Whole Ice Rise
xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];
dx = .5; %[km]
[xx, yy] = meshgrid([xLimits(1):dx:xLimits(2)]*1e3,[yLimits(1):dx:yLimits(2)]*1e3);

% %Promontory
% xLimits = [-3.71e2 -2.95e2];
% yLimits = [-5.32e2 -4.66e2];

%% calc flow
hh = bedmachine_interp('thickness',xx,yy);
% hh = bedmap2_interp(xx,yy,'thickness');

% u_interp = scatteredInterpolant(data2.xy(:,1),data2.xy(:,2),data2.u,'linear','none');
% v_interp = scatteredInterpolant(data2.xy(:,1),data2.xy(:,2),data2.v,'linear','none');
% uu_i = u_interp(xx,yy)*3.154e7;
% vv_i = u_interp(xx,yy)*3.154e7;
% [uui_x,~] = gradient(uu_i.*hh,dx*1e3);
% [~,vvi_y] = gradient(vv_i.*hh,dx*1e3);
% divFlow = uui_x+vvi_y;

divFlow = (data2.A*(data2.u .* data2.h)+ data2.B*(data2.v .* data2.h))*3.154e7;

[uu, vv] = measures_interp('velocity',xx,yy);
uu_b = imgaussfilt(uu.*hh,2/dx);
vv_b = imgaussfilt(vv.*hh,2/dx);
[uu_x,~] = gradient(uu_b,dx*1e3);
[~,vv_y] = gradient(vv_b,dx*1e3);
divFlowMeasures = (uu_x+vv_y);


uu_sm = imgaussfilt(uu,2/dx);
vv_sm = imgaussfilt(vv,2/dx);
hh_sm = imgaussfilt(hh,2/dx);
[uu_xh,~] = gradient(uu_sm,dx*1e3);
[~,vv_yh] = gradient(vv_sm,dx*1e3);
[hh_x,hh_y] = gradient(hh_sm,dx*1e3);

bedComp = uu_sm .* hh_x + vv_sm .* hh_y;
velComp = uu_xh.*hh_sm + vv_yh.*hh_sm;

figure('Position',[300 300 1800 460])
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
ax1 = nexttile(1);
surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - bedComp,'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
% c = colorbar;
% c.Label.String = 'Change in Surface Height [m/yr]';
title('bed div')
colormap(ax1,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
xlabel('Easting [km]')
axis equal

ax2 = nexttile(2);
surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - velComp,'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
% c = colorbar;
% c.Label.String = 'Change in Surface Height [m/yr]';
title('velocity div')
colormap(ax2,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
xlabel('Easting [km]')
axis equal


ax3 = nexttile(3);
surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - divFlowMeasures,'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
% c = colorbar;
% c.Label.String = 'Change in Surface Height [m/yr]';
title('Mass Balannce dh/dt')
colormap(ax3,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
xlabel('Easting [km]')
axis equal






%% Plotting
figure('Position',[300 300 1800 460])
tiledlayout(1,4, 'Padding', 'compact', 'TileSpacing', 'compact');
% ax1 = nexttile(1);
% surf(xx/1e3,yy/1e3,zeros(size(xx)),divFlowMeasures,'edgecolor','none')
% hold on 
% contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
% view(2)
% c = colorbar;
% c.Label.String = '[m/yr]';
% title('Measures')
% colormap(ax1,(cbrewer('div','PiYG',256)));
% caxis(cValues);
% xlim(xLimits)
% ylim(yLimits)
% ylabel('Northing [km]')
% xlabel('Easting [km]')

ax2 = nexttile(1);
% surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - divFlow,'edgecolor','none')
trisurf(data2.t_c,data2.xy_c(:,1)/1e3,data2.xy_c(:,2)/1e3,zeros(size(divFlow)),...
    Acc(data2.xy_c(:,1),data2.xy_c(:,2))*3.154e7 - divFlow,...
    'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
% c = colorbar;
% c.Label.String = 'Change in Surface Height [m/yr]';
title('Modeled dh/dt yr 0')
colormap(ax2,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
ylabel('Northing [km]')
xlabel('Easting [km]')
axis equal


baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt50SpeedUp0.mat";
data2 = load(baseFile);
divFlow = (data2.A*(data2.u .* data2.h)+ data2.B*(data2.v .* data2.h))*3.154e7;


ax2 = nexttile(2);
% surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - divFlow,'edgecolor','none')
trisurf(data2.t_c,data2.xy_c(:,1)/1e3,data2.xy_c(:,2)/1e3,zeros(size(divFlow)),...
    Acc(data2.xy_c(:,1),data2.xy_c(:,2))*3.154e7 - divFlow,...
    'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
% c = colorbar;
% c.Label.String = 'Change in Surface Height [m/yr]';
title('Modeled dh/dt yr 50')
colormap(ax2,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
xlabel('Easting [km]')
axis equal


ax3 = nexttile(3);
surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - imgaussfilt(divFlowMeasures,4),'edgecolor','none')
% surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - divFlowMeasures,'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
% c = colorbar;
% c.Label.String = 'Change in Surface Height [m/yr]';
title('Mass Balannce dh/dt')
colormap(ax3,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
xlabel('Easting [km]')
axis equal

ax4 = nexttile(4);
surf(xx/1e3,yy/1e3,zeros(size(xx)),data2.dhdt_interp(xx,yy),'edgecolor','none')
hold on
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
c = colorbar;
c.Label.String = 'Change in Surface Height [m/yr]';
title('Directly Observed dh/dt')
colormap(ax4,cMapToUse);
caxis(cValues);
xlim(xLimits)
ylim(yLimits)
xlabel('Easting [km]')
axis equal


% figure
% surf(xx/1e3,yy/1e3,zeros(size(xx)),sqrt(uu.^2 + vv.^2),'edgecolor','none')
% view(2)
% colorbar
% c = colorbar;
% c.Label.String = '[m/yr]';
% set(gca,'ColorScale','log')
% title('Obs Speed')

fig = gcf;
setFontSize(24)
labelTiledLayout(fig, 18)

if(saveFigs)
    savePng("figs/paper/" + groupName + fig.Number);
%     saveVect("figs/fig_groupName" + fig.Number);
end

