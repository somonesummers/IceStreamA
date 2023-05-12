function [] = plotThickness(data1,n,ax)
    ftsize = 20;
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
            0*data1.h_s,data1.h...
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
        c.Label.String = 'Ice Thickness [m]';
    end
    title("\Delta" + "z = " + data1.thin_m + "m",'FontSize',ftsize)
    caxis([10^1 10^3.6])     
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    c.FontSize = ftsize;
    view(2)
    axis equal
%     xlim([-4.5e5 -2.5e5])
%     ylim([-5.75e5 -3.5e5])
end

