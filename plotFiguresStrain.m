clear
close all
saveFigs = false;

if(saveFigs)
    disp("Please confirm you'd like to save figures");
    pause()
end

baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
groupName = 'ISSM_N';
cases = [10,0];
speeds = [1,0];

%Whole Ice Rise
xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];

%Promontory
% xLimits = [-3.71e2 -2.95e2];
% yLimits = [-5.32e2 -4.66e2];

file = "data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0";
data1 = load("data/spdChange/" + file + ".mat");

figure('Position',[300 300 1800 733])
tiledlayout(numel(speeds),numel(cases), 'Padding', 'tight', 'TileSpacing', 'tight');


for i = 1:numel(speeds)
    for j = 1:numel(cases)
        if(isfile(strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"DhDt0","DhDt" + cases(j))))
            data1 = load(strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"DhDt0","DhDt" + cases(j)));
            ax = nexttile(j + (i-1)*numel(cases));  
            plotLogStrain(data1,3,ax)
        
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
            xlim(xLimits)
            ylim(yLimits)
        else
            warning("File not found: " + strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"Thin0","Thin" + cases(j)));
        end
    end
end