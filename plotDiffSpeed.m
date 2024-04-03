function [] = plotDiffSpeed(data1, data2, n,ax)
%plots speed difference between data1 (exp) and data2 (control)
    ftsize = 20;
    trisurf(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            0*data1.h_s,(sqrt(data1.u.^2 + data1.v.^2)*3.154E7...
            - sqrt(data2.u.^2 + data2.v.^2)*3.154E7)...
               ,'edgecolor','none','facecolor','interp');
    hold on
    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H1] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
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
        [~, H2] = tricontour(data2.t,data2.xy(:,1)/1e3,data2.xy(:,2)/1e3,...
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
        ylabel('Northing [km]')
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
%     xlim([-4.5e5 -2.5e5])
%     ylim([-5.75e5 -3.5e5])
end

