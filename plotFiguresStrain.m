clear
close all
saveFigs = false;

if(saveFigs)
    disp("Saving figures, make sure you're not overwriting anything you like");
%     pause()
end

baseFile = "data/spdChange/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
groupName = 'ISSM_linear';
cases = [50:-20:0,0];
speeds = [5:-2:0,0];

%Whole Ice Rise
xLimits = [-5.1e2 -2.5e2];
yLimits = [-5.75e2 -3.5e2];

%Promontory
% xLimits = [-3.71e2 -2.95e2];
% yLimits = [-5.32e2 -4.66e2];

file = "data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0";
data1 = load("data/spdChange/" + file + ".mat");





for i = 1:3
    figure('Position',[300 300 1800 500])
    tiledlayout(numel(1),numel(cases), 'Padding', 'tight', 'TileSpacing', 'tight');
    
    strainType = i; % 1 lat, 2 trans, 3 shear

    if(strainType == 1)
        headTitle = "Lon";
    elseif(strainType == 2)
        headTitle = "Trans";
    elseif(strainType == 3)
        headTitle = "Shear";
    else
        error("Bad Strain Type")
    end
    for j = 1:numel(cases)
        if(isfile(strrep(strrep(baseFile,"Up0","Up" + speeds(j)),"DhDt0","DhDt" + cases(j))))
%             tic
            data1 = load(strrep(strrep(baseFile,"Up0","Up" + speeds(j)),"DhDt0","DhDt" + cases(j)),...
                'xy','u','v','h','pv','t');
            ax = nexttile(j);  
            plotLogStrain(data1,strainType,ax)
        
            if(i ==1)
                title("\Delta z = " + cases(j) + " " + (10+speeds(j))*10 + "%")
            end
            if(j == numel(cases)) 
                c = colorbar;
                    c.Label.String = 'Strain [yr^{-1}]';
    %                 c.FontSize = 18;
            end
            xlabel("")
            ylabel("")
            xlim(xLimits)
            ylim(yLimits)
%             toc
        else
            warning("File not found: " + strrep(strrep(baseFile,"Up0","Up" + speeds(i)),"Thin0","Thin" + cases(j)));
        end
        setFontSize(18)
    end


    sgtitle(headTitle,'Fontsize',24)

    if(saveFigs)
        savePng("figs/paper/Strain" + headTitle +"_" + groupName);
        close
    end

end