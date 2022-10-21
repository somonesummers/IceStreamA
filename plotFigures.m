%% Plot all figures
clear
close all
saveFigs = false;

%% Cases of thickness

cases = [-50,-20,0,20,50];
figure('Position',[300 300 1000 1000])
tiledlayout(2,numel(cases), 'Padding', 'none', 'TileSpacing', 'tight');

data2 = load("data/data_gridRefinedXSM025ISSM Shift5bedmap0.mat");
% [uu,vv] = measures_interp('velocity',data2.xy(:,1),data2.xy(:,2));
% data2.u = uu/3.154E7;
% data2.v = vv/3.154E7;
for j = 1:numel(cases)
    if(isfile("data/data_gridRefinedXSM025ISSM Shift5bedmap"+cases(j)+".mat"))
        data1 = load("data/data_gridRefinedXSM025ISSM Shift5bedmap"+cases(j)+".mat");    
        ax1 = nexttile(j);  
        plotSpeed(data1,0,ax1);
        if(j == 1)
            ylabel("Northing [m]",'fontsize',18);
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Ice Speed [m/yr]';
            f = gca;
            f.ColorScale = 'log';
            c.FontSize = 18;
        end
        xlabel("");
        title("Case " + cases(j))
        ax2 = nexttile(j+numel(cases));
        plotDiff(data1,data2,0,ax2);
        if(j == 1)
            ylabel("Northing [m]",'fontsize',18);
        end
        if(j == numel(cases))
          	c = colorbar;
            c.Label.String = 'Speed Diff [m/yr]';
            c.FontSize = 18;
        end
        xlabel("Easting [m]",'fontsize',18);
        caxis([-100,100])
        xlabel("");
%         ax3 = nexttile(j+8);
%         plotTau(data1,j,ax3);
    else
        warning("File not found: data/data_ThwaitesBasinGrid" + cases(j) + "0_Uniform_ConstantFlux.mat");
    end
end

if(saveFigs)
    fig = gcf;
    savePng("figs/fig_" + fig.Number);
    saveVect("figs/fig_" + fig.Number);
end


%% Plot xi
ftsize = 20;
load('gridSipleXXSmall3000.mat');
dx = 1e3;
figure('Position',[300 300 700 850])
x = [-5e5:dx:-2e5];
y = [-6e5:dx: -3e5];
[xx,yy] = meshgrid(x,y);
sf_raw =  bedmachine_interp('surface',xx,yy);
sf = imgaussfilt(sf_raw,5e3/dx);
h  =  bedmachine_interp('thickness',xx,yy);
spd        = measures_interp('speed',xx,yy);
[sx ,  sy] = gradient(sf,dx,dx);

% p = surf(xx,yy,zeros(size(h)),h./(sqrt(sx.^2+sy.^2)*200e3));
p = surf(xx,yy,zeros(size(h)),spd);
hold on
plot(pv(:,1),pv(:,2),'k-','LineWidth', 4)
contour(x,y,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(x,y,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(x,y,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(x,y,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
% [C,hh] = contour(x,y,h./(sqrt(sx.^2+sy.^2)*200e3),[.1,.3,1,3,10], 'r-','HandleVisibility','off');  
% clabel(C,hh)
% title('\xi factor')
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
% colormap(cmocean('curl'))
set(p, 'edgecolor', 'none');
caxis([10^1 10^3.6]) 
axis equal
ylabel('Northing [m]')
xlabel('Easting [m]')
c = colorbar;
c.Label.String = 'Ice Speed [m/yr]';
c.FontSize = ftsize;   
view(2)
f = gca;
f.XAxis.FontSize = ftsize-2;
f.YAxis.FontSize = ftsize-2;
view(2)
axis equal
xlim([-5e5, -2e5])
ylim([-6e5, -3e5])
if(saveFigs)
    fig = gcf;
    savePng("figs/fig_" + fig.Number);
    saveVect("figs/fig_" + fig.Number);
end


%% Functions
function [] = plotDiff(data1, data2, n,ax)
    ftsize = 20;
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
            0*data1.h_s,(sqrt(data1.u.^2 + data1.v.^2)*3.154E7...
            - sqrt(data2.u.^2 + data2.v.^2)*3.154E7)...
               ,'edgecolor','none','facecolor','interp');
    hold on
    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H1] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
            sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H1)
        H1(j).EdgeColor = rgb('black');
        H1(j).LineStyle = '--';
        H1(j).LineWidth = 3;
        end
    else
        warning('Speed 1 seems wrong, no 30 m/a contour')
    end

    if(max(sqrt(data2.u.^2 + data2.v.^2)*3.154E7) > 30)
        [~, H2] = tricontour(data2.t,data2.xy(:,1),data2.xy(:,2),...
            sqrt(data2.u.^2 + data2.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H2)
        H2(j).EdgeColor = rgb('grey');
        H2(j).LineStyle = '--';
        H2(j).LineWidth = 3;
        end
    else
        warning('Speed 2 seems wrong, no 30 m/a contour')
    end
    
    
    
%     CT = flipud(cbrewer('div','RdBu',256));
    CT = cmocean('balance');
    colormap(ax, CT)
%     colormap(ax, 'redblue')
    if(n == 2 || n < 0)
        ylabel('Northing [m]')
    end
    xlabel('Easting [m]')
    if(n == 4 || n < 0)
        c = colorbar;
        c.Label.String = '\Delta Ice Speed [m/yr]';
        c.FontSize = ftsize;
    end
%     title(data1.thin_m)
    f = gca;
%     caxis([-750 750]) %observed diffs red blue
    caxis([-325 325]) %observed diffs red blue
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    view(2)
    axis equal
    xlim([-4.5e5 -2.5e5])
    ylim([-5.75e5 -3.5e5])
end

function [] = plotSpeed(data1,n,ax)
    ftsize = 20;
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
            0*data1.h_s,sqrt(data1.u.^2 + data1.v.^2)*3.154E7...
               ,'edgecolor','none','facecolor','interp');
    hold on
    
    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
            sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H)
            if(n == 1 || n == -1)
                H(j).EdgeColor = rgb('gray');
            else
                H(j).EdgeColor = rgb('black');
            end
            H(j).LineStyle = '--';
            H(j).LineWidth = 3;
        end
    else
        warning('Speed 1 seems wrong, no 30 m/a contour')
    end

    if(n == 1 || n < 0)
        ylabel('Northing [m]')
    end
    xlabel('Easting [m]')
    if(n == 4 || n < 0)
        c = colorbar;
        c.Label.String = 'Ice Speed [m/yr]';
    end
    title("\Delta" + "z = " + data1.thin_m + "m",'FontSize',ftsize)
    caxis([10^1 10^3.6])     
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    f.ColorScale = 'log';
    c.FontSize = ftsize;
    view(2)
    axis equal
    xlim([-4.5e5 -2.5e5])
    ylim([-5.75e5 -3.5e5])
end

function [] = plotThinning(data1,data2,n,ax)
    load('thinningrates.mat','thinInterpConic');
    ftsize = 20;
    if (n == 1)
        trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
            0*data1.h,(-data1.thin_m*ones(size(data1.h)))...
               ,'edgecolor','none','facecolor','interp');
    elseif (n == 2)
        trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
            0*data1.h,data1.thin_m*thinInterpConic(data1.xy(:,1),data1.xy(:,2))...
               ,'edgecolor','none','facecolor','interp');
    end
    hold on
    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
            sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H)
        H(j).EdgeColor = rgb('black');
        H(j).LineStyle = '--';
        H(j).LineWidth = 3;
        end
    else
        warning('Speed 1 seems wrong, no 30 m/a contour')
    end

    if(max(sqrt(data2.u.^2 + data2.v.^2)*3.154E7) > 30)
        [~, H] = tricontour(data2.t,data2.xy(:,1),data2.xy(:,2),...
            sqrt(data2.u.^2 + data2.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H)
        H(j).EdgeColor = rgb('grey');
        H(j).LineStyle = '--';
        H(j).LineWidth = 3;
        end
    else
        warning('Speed 2 seems wrong, no 30 m/a contour')
    end
    CT2 = flipud(cbrewer('seq','Oranges',14));
    colormap(ax, CT2)
    ylabel('Northing [m]')
    xlabel('Easting [m]')
    c = colorbar;
    c.Label.String = 'Applied \Delta H [m]';
    c.FontSize = ftsize;   
    caxis([-70 0])
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    view(2)
    axis equal
%     xlim([-1.56e6 -1.38e6])
%     ylim([-5.75e5 -3e5])
end

function [] = plotFlux(data1, data2, ax)
    plot(data1.h_init(data1.xy(data1.dwnSt_bound,1),data1.xy(data1.dwnSt_bound,2))...
        .*sqrt(data1.u(data1.dwnSt_bound).^2+data1.v(data1.dwnSt_bound).^2)*3.154E7...
        ,'k')
    hold on
    plot(data2.h_init(data1.xy(data1.dwnSt_bound,1),data1.xy(data1.dwnSt_bound,2))...
        .*sqrt(data2.u(data1.dwnSt_bound).^2+data2.v(data1.dwnSt_bound).^2)*3.154E7...
        ,'--','color',rgb('gray'),'linewidth',2)
    
    plot(data1.h_init(data1.xy(data1.upSt_bound,1),data1.xy(data1.upSt_bound,2))...
        .*sqrt(data1.u(data1.upSt_bound).^2+data1.v(data1.upSt_bound).^2)*3.154E7...
        ,'k')
    plot(data2.h_init(data1.xy(data1.upSt_bound,1),data1.xy(data1.upSt_bound,2))...
        .*sqrt(data2.u(data1.upSt_bound).^2+data2.v(data1.upSt_bound).^2)*3.154E7...
        ,'--','color',rgb('gray'),'linewidth',2)
    
end

function[] = plotTau(data1,n,ax)
    ftsize = 20;
    load Dawn.mat
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),data1.tau_c(data1.xy(:,1),...
        data1.xy(:,2),data1.u,data1.v,data1.grounded)./norms([data1.u,data1.v],2,2)/1e3,...
       'edgecolor','none','facecolor','interp')
%     hold on
%     trisurf(t,xy(:,1),xy(:,2),-10e2*ones(size(xy(:,1))),...
%            'edgecolor','black','facecolor','none')
    caxis([0 250]);
    colormap(ax, (Cmap/255.0))
%     title("Basal Strength",'FontSize',ftsize+2)
    if(n == 1 || n <= 0)
        ylabel('Northing [m]')
    end
    xlabel('Easting [m]')
    if(n == 4 || n <= 0)
        c = colorbar;
        c.Label.String = 'Basal Strength [kPa]';
    end
    view(2)
    axis equal
    % axis off
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    c.FontSize = ftsize;
%     xlim([-1.56e6 -1.38e6])
%     ylim([-5.75e5 -3e5])
end

function [] = plotDiffStress(data1, data2, ax)
    ftsize = 20;
    trisurf(data1.t_c,data1.xy_c(:,1),data1.xy_c(:,2),...
            0*data1.h_av,100*(data1.df-data2.df)./data2.df...
               ,'edgecolor','none','facecolor','interp');
    hold on

    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
            sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H)
        H(j).EdgeColor = rgb('black');
        H(j).LineStyle = '--';
        H(j).LineWidth = 3;
        end
    else
        warning('Speed 1 seems wrong, no 30 m/a contour')
    end

    if(max(sqrt(data2.u.^2 + data2.v.^2)*3.154E7) > 30)
        [~, H] = tricontour(data2.t,data2.xy(:,1),data2.xy(:,2),...
            sqrt(data2.u.^2 + data2.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H)
        H(j).EdgeColor = rgb('grey');
        H(j).LineStyle = '--';
        H(j).LineWidth = 3;
        end
    else
        warning('Speed 2 seems wrong, no 30 m/a contour')
    end

    CT2 = flipud(cbrewer('div','PuOr',256));
%     title(data1.thin_m)
    colormap(ax, CT2);
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    c = colorbar;
    caxis([-10 10]) %observed diffs red blue
    c.Label.String = '\Delta Driving Stess Change [%]';
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    c.FontSize = ftsize;
    view(2)
    axis equal
%     xlim([-1.56e6 -1.38e6])
%     ylim([-5.75e5 -3e5])
   
end