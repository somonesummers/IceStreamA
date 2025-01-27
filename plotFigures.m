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
groupName = 'massBlanceDhDt';
cases = [30,20,10,0]; %thinning cases
titleCases = ["-9m elevation","-6m elevation","-3m elevation","-0m elevation"];
% cases = [5:-1:-1]; %speed cases
figure('Position',[300 300 1200 733])
tiledlayout(3,numel(cases), 'Padding', 'compact', 'TileSpacing', 'compact');

% baseFile = "data/DhDt/data_NgridFlowRiseA02ISSMDhDt0case.mat";
% baseFile = "data/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
% baseFile = "data/spdChange/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
baseFile = "data/datagridFlowRiseA02UpBCISSMMb_DhDt0SpeedUp0.mat";


data2 = load(baseFile);

%Whole Ice Rise
xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];

% %Promontory
% xLimits = [-3.71e2 -2.95e2];
% yLimits = [-5.32e2 -4.66e2];



for j = 1:numel(cases)
%     baseFile = "data/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
    newFile = strrep(baseFile,"Dt0","Dt" + cases(j));
%     newFile = strrep(baseFile,"Up0","Up" + cases(j));
%     baseFile = strrep(newFile,"_N","");
    data2 = load(baseFile);
    if(isfile(newFile))
        data1 = load(newFile);    
        ax1 = nexttile(j);  
        plotSpeed(data1,0,ax1);
%         plotNef(data1,0,ax1);
        axis equal
        xlim(xLimits)
        ylim(yLimits)
        xticklabels([])
        if(j == 1)
            ylabel("Northing [km]",'fontsize',18);
        else
            yticklabels([])
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Speed [m/yr]';
            c.FontSize = 18;
        end
        xlabel("");
        title(titleCases(j))
%         title((cases(j)+10)*10 +"% Velocity")
        if(j == 3)
%             title(groupName)
        end
        ax2 = nexttile(j+numel(cases));
%         plotThickness(data1,0,ax2);
%         plotDiffHeight(data1,data2,0,ax2);
        plotDiffSpeed(data1,data2,0,ax2);
        xlim(xLimits)
        ylim(yLimits)
        xticklabels([])
%         plotTau(data1,0,ax2);
        if(j == 1)
            ylabel("Northing [km]",'fontsize',18);
        else
            yticklabels([])
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Speed Diff [m/yr]';
            c.FontSize = 18;
        end
        xlabel("")
        xticklabels([])
        
        ax3 = nexttile(j+2*numel(cases));
        plotStress(data1,0,ax3)
        xlim(xLimits)
        ylim(yLimits)

        title("")
%         plotDiffHeight(data1,data2,0,ax3);
%         plotTau(data1,0,ax3);
%         plotDiffTau(data1,data2,0,ax3);
        if(j == 1)
            ylabel("Northing [km]",'fontsize',18);
        else
            yticklabels([])
        end
        if(j == numel(cases))
          	c = colorbar;
%             c.Label.String = 'Surf Height Diff [m]';
%             c.Label.String = 'Strength [kPa]';
            c.Label.String = 'Surface Stress [Pa]';
            c.FontSize = 18;
        end
    else
        warning("File not found:" + cases(j));
    end
end

if(saveFigs)
    fig = gcf;
    labelTiledLayout(fig, 18)
    savePng("figs/paper/" + groupName + fig.Number);
%     saveVect("figs/fig_groupName" + fig.Number);
end

