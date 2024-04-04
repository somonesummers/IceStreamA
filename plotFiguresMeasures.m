%% Plot all figures for 1D change, This is best plotting option
clear
% close all
addpath lib/
saveFigs = true;

if(saveFigs)
    disp("Please confirm you'd like to save figures");
    pause()
end
%% Cases of thickness
% time = load('Golledge21_GRL_T1_thick_22mar23_v2_Paul.mat','time');
groupName = 'ISSM_Base';
figure('Position',[300 300 1800 350])
tiledlayout(1,4, 'Padding', 'tight', 'TileSpacing', 'tight');

% baseFile = "data/DhDt/data_NgridFlowRiseA02ISSMDhDt0case.mat";
% baseFile = "data/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";

data2 = load(baseFile);

[uu,vv] = measures_interp('velocity',data2.xy(:,1),data2.xy(:,2));
data2.u = uu/3.154E7;
data2.v = vv/3.154E7;

%Whole Ice Rise
xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];

% %Promontory
% xLimits = [-3.71e2 -2.95e2];
% yLimits = [-5.32e2 -4.66e2];



newFile = baseFile;

if(isfile(newFile))
    data1 = load(newFile);    
    
    ax1 = nexttile(1);  
    plotSpeed(data2,0,ax1);
    axis equal
    xlim(xLimits)
    ylim(yLimits)
    
    ylabel("Northing [km]",'fontsize',18);
    xlabel("Easting [km]",'fontsize',18);
    c = colorbar;
    c.Label.String = 'Speed [m/yr]';
    c.FontSize = 18;
    title('MEaSUREs Speed','fontsize',20)
    
    
    ax2 = nexttile(2);  
    plotSpeed(data1,0,ax2);
    axis equal
    xlim(xLimits)
    ylim(yLimits)
    
    ylabel("Northing [km]",'fontsize',18);
    xlabel("Easting [km]",'fontsize',18);
    c = colorbar;
    c.Label.String = 'Speed [m/yr]';
    c.FontSize = 18;
    title('Model Speed','fontsize',20)
    
    ax3 = nexttile(3);
    plotDiffSpeed(data1,data2,0,ax3);
    xlim(xLimits)
    ylim(yLimits)
    ylabel("Northing [km]",'fontsize',18);
    xlabel("Easting [km]",'fontsize',18);
    c = colorbar;
    c.Label.String = 'Speed Diff [m/yr]';
    c.FontSize = 18;
    title('Model Misfit','fontsize',20)

    ax4 = nexttile(4);
    plotTau(data1,0,ax4)
    xlim(xLimits)
    ylim(yLimits)
    ylabel("Northing [km]",'fontsize',18);
    xlabel("Easting [km]",'fontsize',18);
    c = colorbar;
    c.Label.String = 'Basal Strength [kPa]';
    c.FontSize = 18;
    title('Model \tau','fontsize',20)
else
    warning("File not found: " + newFile);
end

if(saveFigs)
    fig = gcf;
    labelTiledLayout(fig, 18)
    savePng("figs/paper/" + groupName + fig.Number);
%     saveVect("figs/fig_groupName" + fig.Number);
end

