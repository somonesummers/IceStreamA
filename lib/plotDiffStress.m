
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

