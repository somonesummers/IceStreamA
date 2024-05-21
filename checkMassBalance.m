%% Plot all figures for 1D change, This is best plotting option
clear
close all
addpath lib/
saveFigs = false;

if(saveFigs)
    disp("Please confirm you'd like to save figures");
    pause()
end

[Acc, T_s] = loadALBMAP();


%% Cases of thickness
groupName = 'ISSM_N_thinning';
% cases = [5:-1:-1]; %speed cases
figure('Position',[300 300 1800 733])
tiledlayout(1,4, 'Padding', 'compact', 'TileSpacing', 'compact');
baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";

data2 = load(baseFile);

%Whole Ice Rise
xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];
dx = .5; %[km]
[xx, yy] = meshgrid([xLimits(1):dx:xLimits(2)]*1e3,[yLimits(1):dx:yLimits(2)]*1e3);

% %Promontory
% xLimits = [-3.71e2 -2.95e2];
% yLimits = [-5.32e2 -4.66e2];

divFlow = (data2.A*(data2.u .* data2.h)+ data2.B*(data2.v .* data2.h))*3.154e7;
h = bedmachine_interp('thickness',xx,yy);
[uu, vv] = measures_interp('velocity',xx,yy);
uu_b = imgaussfilt(uu.*h,2/dx);
vv_b = imgaussfilt(vv.*h,2/dx);
[uu_x,~] = gradient(uu_b,dx*1e3);
[~,vv_y] = gradient(vv_b,dx*1e3);
divFlowMeasures = (uu_x+vv_y);

ax1 = nexttile(1);
surf(xx/1e3,yy/1e3,zeros(size(xx)),divFlowMeasures,'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
c = colorbar;
c.Label.String = '[m/yr]';
title('Measures')
colormap(ax1,flipud(cbrewer('div','RdBu',256)));
caxis([-2 2]);
xlim(xLimits)
ylim(yLimits)

ax2 = nexttile(2);
trisurf(data2.t_c,data2.xy_c(:,1)/1e3,data2.xy_c(:,2)/1e3,zeros(size(divFlow)),...
    Acc(data2.xy_c(:,1),data2.xy_c(:,2))*3.154e7 - divFlow,...
    'edgecolor','none','facecolor','interp')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
c = colorbar;
c.Label.String = '[m/yr]';
title('Atm - Model Div')
colormap(ax2,flipud(cbrewer('div','RdBu',256)));
caxis([-.7 .7]);
xlim(xLimits)
ylim(yLimits)

ax3 = nexttile(3);
surf(xx/1e3,yy/1e3,zeros(size(xx)),Acc(xx,yy)*3.154e7 - divFlowMeasures,'edgecolor','none')
hold on 
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
c = colorbar;
c.Label.String = '[m/yr]';
title('Atm Model - Measures Div')
colormap(ax3,flipud(cbrewer('div','RdBu',256)));
caxis([-.7 .7]);
xlim(xLimits)
ylim(yLimits)

ax4 = nexttile(4);
surf(xx/1e3,yy/1e3,zeros(size(xx)),data2.dhdt_interp(xx,yy),'edgecolor','none')
hold on
contour(xx/1e3,yy/1e3,sqrt(uu.^2 + vv.^2),[30 30],'k')
view(2)
c = colorbar;
c.Label.String = '[m/yr]';
title('Obs changes')
colormap(ax4,flipud(cbrewer('div','RdBu',256)));
caxis([-.7 .7]);
xlim(xLimits)
ylim(yLimits)



% figure
% surf(xx/1e3,yy/1e3,zeros(size(xx)),sqrt(uu.^2 + vv.^2),'edgecolor','none')
% view(2)
% colorbar
% c = colorbar;
% c.Label.String = '[m/yr]';
% set(gca,'ColorScale','log')
% title('Obs Speed')


if(saveFigs)
    fig = gcf;
    labelTiledLayout(fig, 18)
    savePng("figs/paper/" + groupName + fig.Number);
%     saveVect("figs/fig_groupName" + fig.Number);
end

