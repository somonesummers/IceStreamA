function[] = plotNef(data1,n,ax)
    ftsize = 20;
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),zeros(size(data1.N_ef)),data1.N_ef,...
       'edgecolor','none','facecolor','interp')
    if(n == 1 || n < 0)
        ylabel('Northing [m]')
    end
    xlabel('Easting [m]')
    if(n == 4 || n < 0)
        c = colorbar;
        c.Label.String = 'Effect Pressure [Pa]';
    end
    hold on
    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H2] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
            sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
        for j = 1:numel(H2)
        H2(j).EdgeColor = rgb('black');
        H2(j).LineStyle = '-';
        H2(j).LineWidth = 1;
        end
    else
        warning('Speed 2 seems wrong, no 30 m/a contour')
    end
    view(2)
    axis equal
    CT = cmocean('balance');
    caxis([.4 1.6])
    colormap(ax, CT)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    c.FontSize = ftsize;
%     xlim([-1.56e6 -1.38e6])
%     ylim([-5.75e5 -3e5])
end

