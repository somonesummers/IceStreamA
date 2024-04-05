%% Plot all figures for 1D change, This is best plotting option
clear
% close all
addpath lib/
saveFigs = false;

%% Cases of thickness
groupName = 'ISSM';
cases = [50,10,0];
figure('Position',[300 300 1300 680])
tiledlayout(3,numel(cases), 'Padding', 'none', 'TileSpacing', 'tight');

baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp2.mat";

data2 = load(baseFile);

for j = 1:numel(cases)
    newFile = strrep(baseFile,"Dt0","Dt" + cases(j));
    data2 = load(baseFile);
    
    if(isfile(newFile))
        data1 = load(newFile);    
        ax1 = nexttile(j);  
        plotStress(data1,0,ax1)

        if(j == 1)
            ylabel("Northing [m]",'fontsize',18);
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Effective Stress [Pa]';
            c.FontSize = 18;
        end
        xlabel("");
        title(cases(j) + " years ago")
%         title("Year " + time.time(cases(j)))
        if(j == 3)
%             title(groupName)
        end
        ax2 = nexttile(j+numel(cases));
        
        plotStressComp(data1,1,0,ax2);
        
        if(j == 1)
            ylabel("Northing [m]",'fontsize',18);
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Stress 1 [Pa]';
            c.FontSize = 18;
        end
        xlabel("Easting [m]",'fontsize',18);
        xlabel("");
        
        ax3 = nexttile(j+2*numel(cases));
        plotStressComp(data1,2,0,ax3);
        
        if(j == 1)
            ylabel("Northing [m]",'fontsize',18);
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Stress 2 [Pa]';
            c.FontSize = 18;
        end
    else
        warning("File not found:" + cases(j));
    end
end

if(saveFigs)
    fig = gcf;
    savePng("figs/fig_" + groupName + fig.Number);
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