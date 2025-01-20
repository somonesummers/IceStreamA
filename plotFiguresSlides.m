%% Plot all figures
clear
% close all
addpath lib/
saveFigs = false;

%% Cases of thickness
% time = load('Golledge21_GRL_T1_thick_22mar23_v2_Paul.mat','time');
groupName = 'SlidesNoHAF';
cases = [40,0];
figure('Position',[300 300 1000 680])
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'tight');

baseFile = "data/DhDt/data_NgridFlowRiseA02ISSMDhDt0Speedup0.mat";
% baseFile = "data/NoHAF/datagridFlowRiseA02ISSMDhDt0SpeedUp0.mat";
% baseFile = "data/data_NgridFlowRiseA05ISSMGoll61case.mat";

% Below utilizes sshfs to directly use files on server. It is slow, but
% allows me to not use up space on the laptop hardrive.
% baseFile = "~/sherlock_home/IceStreamA/data/data_NgridFlowRiseA02ISSMDhDt0Tau1_30.mat";

% Below is from external harddrive, faster, also not this HD
% baseFile = '/Volumes/Extreme SSD/IceStreamAData/datagridFlowRiseA02ISSMDhDt0SpeedUp0.mat';
% baseFile = '/Volumes/Extreme SSD/IceStreamAData/data_NgridFlowRiseA02ISSMDhDt0case.mat';
% baseFile = '/Volumes/Extreme SSD/IceStreamAData/data_NgridFlowRiseA02ISSMGoll21caseMay15.mat';
% baseFile = '/Volumes/Extreme SSD/IceStreamAData/data_NgridFlowRiseA02ISSMGoll60case.mat';
data2 = load(baseFile);
% [uu,vv] = measures_interp('velocity',data2.xy(:,1),data2.xy(:,2));z`
% data2.u = uu/3.154E7;
% data2.v = vv/3.154E7;
newFile = strrep(baseFile,"Dt0","Dt" + cases(1));
data1 = load(newFile);    
ax1 = nexttile(1);  
plotSpeed(data1,0,ax1);
ylabel("Northing [m]",'fontsize',18);
title(cases(1)*.5 + " meters case")

ax1 = nexttile(2);  
plotSpeed(data2,0,ax1);
c = colorbar;
c.Label.String = 'Speed [m/yr]';
c.FontSize = 18;
title(cases(2)*.5 + " meters case")

ax2 = nexttile(3);
%         plotThickness(data1,0,ax2);
plotDiffSpeed(data1,data2,0,ax2);
%         plotTau(data1,0,ax2);

ylabel("Northing [m]",'fontsize',18);
c = colorbar;
c.Label.String = 'Speed Diff [m/yr]';
c.FontSize = 18;
xlabel("Easting [m]",'fontsize',18);


ax3 = nexttile(4);
plotDiffHeight(data1,data2,0,ax3);
ylabel("Northing [m]",'fontsize',18);

c = colorbar;
c.Label.String = 'Surf Height Diff [m]';
%             c.Label.String = 'Strength [kPa]';
c.FontSize = 18;

if(saveFigs)
    fig = gcf;
    savePng("figs/fig_" + groupName);
%     saveVect("figs/fig_groupName" + fig.Number);
end


%% Plot speed
% ftsize = 20;
% load('grids/gridFlowRiseA02.mat');
% dx = 1e3;
% figure('Position',[300 300 700 850])
% x = [-5e5:dx:-2e5];
% y = [-6e5:dx: -3e5];
% [xx,yy] = meshgrid(x,y);
% sf_raw =  bedmachine_interp('surface',xx,yy);
% sf = imgaussfilt(sf_raw,5e3/dx);
% h  =  bedmachine_interp('thickness',xx,yy);
% spd        = measures_interp('speed',xx,yy);
% [sx ,  sy] = gradient(sf,dx,dx);
% 
% % p = surf(xx,yy,zeros(size(h)),h./(sqrt(sx.^2+sy.^2)*200e3));
% p = surf(xx,yy,zeros(size(h)),spd);
% hold on
% plot(pv(:,1),pv(:,2),'k-','LineWidth', 4)
% contour(x,y,spd, [10, 10] , 'k:','HandleVisibility','off')
% contour(x,y,spd, [30, 30] , 'k--','HandleVisibility','off')
% contour(x,y,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(x,y,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
% % [C,hh] = contour(x,y,h./(sqrt(sx.^2+sy.^2)*200e3),[.1,.3,1,3,10], 'r-','HandleVisibility','off');  
% % clabel(C,hh)
% % title('\xi factor')
% f = gca;
% f.ColorScale = 'log';
% view(2)
% colorbar
% % colormap(cmocean('curl'))
% set(p, 'edgecolor', 'none');
% caxis([10^1 10^3.6]) 
% axis equal
% ylabel('Northing [m]')
% xlabel('Easting [m]')
% c = colorbar;
% c.Label.String = 'Ice Speed [m/yr]';
% c.FontSize = ftsize;   
% view(2)
% f = gca;
% f.XAxis.FontSize = ftsize-2;
% f.YAxis.FontSize = ftsize-2;
% view(2)
% axis equal
% xlim([-5e5, -2e5])
% ylim([-6e5, -3e5])
% if(saveFigs)
%     fig = gcf;
%     savePng("figs/fig_" + groupName + fig.Number);
% %     saveVect("figs/fig_" + fig.Number);
% end