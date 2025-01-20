%% Plot figure for 2D change (ie BC speed and thickness)
clear
close all
saveFigs = false;
addpath('data');

if(saveFigs)
    disp("Please confirm you'd like to save figures"); %this has 
    pause()
end

%% Cases of Thickness/Speed
groupName = 'ISSM_N';
cases = [50:-10:-50];
speeds = [5:-1:-5];
% cases = [10,0];
% speeds = [1,0];


baseFile = "data/spdChange/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
data2 = load(baseFile);
% [uu,vv] = measures_interp('velocity',data2.xy(:,1),data2.xy(:,2));
% data2.u = uu/3.154E7;
% data2.v = vv/3.154E7;

f1 = figure('Position',[300 300 1300 1300]);
tiledlayout(numel(speeds),numel(cases), 'Padding', 'tight', 'TileSpacing', 'none');

for i = 1:numel(speeds)
    for j = 1:numel(cases)

        if(isfile(strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"DhDt0","DhDt" + cases(j))))
            data1 = load(strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"DhDt0","DhDt" + cases(j)));    
            ax1 = nexttile(j + (i-1)*numel(cases));  
            plotDiffSpeed(data1,data2,0,ax1);
            clear ax1
            if(i ==1)
                title("\Delta z = " + cases(j))
            end
            if(j == 1)
                ylabel((10+speeds(i))*10 + "%",'fontsize',18);
            else
                yticklabels([]);
            end
            if(j == numel(cases) && i == round(numel(speeds)/2)) % only plot on center
                c = colorbar;
%                 c.Label.String = 'Speed Diff [m/yr]';
%                 c.FontSize = 18;
            end
            if(i == numel(speeds))
                xlabel("");
            else
                xlabel("")
                xticklabels([]);
            end
            caxis([-100,100])
            xlim([-5.1e2 -2.5e2])
            ylim([-5.75e2 -3.5e2])
    %         ax3 = nexttile(j+8);
    %         plotTau(data1,j,ax3);
        else
            warning("File not found: " + strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"Thin0","Thin" + cases(j)));
        end
    end
end
setFontSize(18)
% labelTiledLayout(f1,20)

f2 = figure('Position',[300 300 1300 1300]);
tiledlayout(numel(speeds),numel(cases), 'Padding', 'tight', 'TileSpacing', 'none');

for i = 1:numel(speeds)
    for j = 1:numel(cases)

        if(isfile(strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"DhDt0","DhDt" + cases(j))))
            data1 = load(strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"DhDt0","DhDt" + cases(j)));    
            ax1 = nexttile(j + (i-1)*numel(cases));  
%             plotStressComp(data1,2,0,ax1);
            plotStress(data1,0,ax1)
            clear ax1 
            title("");
            if(i ==1)
                title("\Delta z = " + cases(j))
            end
            if(j == 1)
                 ylabel((10+speeds(i))*10 + "%",'fontsize',18);
            else
                yticklabels([]);
            end
            if(j == numel(cases) && i == round(numel(speeds)/2)) %only plot on center
                c = colorbar;
%                 c.Label.String = 'Surface Stress [Pa]';
%                 c.FontSize = 18;
            end
            if(i == numel(speeds))
                xlabel("Easting [km]",'fontsize',18);
            else
                xticklabels([]);
                xlabel("");
            end
            xlim([-5.1e2 -2.5e2])
            ylim([-5.75e2 -3.5e2])
%             caxis([-100,100])
    %         ax3 = nexttile(j+8);
    %         plotTau(data1,j,ax3);
        else
            warning("File not found: " + strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"Thin0","Thin" + cases(j)));
        end
    end
end
setFontSize(18)
% labelTiledLayout(f2,20)

if(saveFigs)
    figure(f1);
    savePng("figs/paper/2DSpeed_" + groupName + f1.Number);
    figure(f2);
    savePng("figs/paper/2DStress_" + groupName + f2.Number);
%     saveVect("figs/fig_groupName" + fig.Number);
end

