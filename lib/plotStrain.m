function [] = plotStrain(data1,n,ax)
    ftsize = 20;
    trisurf(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            data1.ep_dot(data1.xy(:,1),data1.xy(:,2))*3.154E7...
               ,'edgecolor','none','facecolor','interp');
    hold on
    
    e_threshold = 5e-3;
    if(max(data1.ep_dot(data1.xy(:,1),data1.xy(:,2))*3.154E7) > e_threshold)
        [~, H] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            data1.ep_dot(data1.xy(:,1),data1.xy(:,2))*3.154E7,[e_threshold e_threshold]);
        for j = 1:numel(H) %there can be multiple contours, so plot them all
            H(j).EdgeColor = rgb('lime green');
            H(j).LineStyle = ':';
            H(j).LineWidth = 2;
        end
    else
        warning('strain seems wrong, no strain threshold contour')
    end
    
%     if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
%         [~, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
%             sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
%         for j = 1:numel(H)
%             if(n == 1 || n == -1)
%                 H(j).EdgeColor = rgb('gray');
%             else
%                 H(j).EdgeColor = rgb('black');
%             end
%             H(j).LineStyle = '--';
%             H(j).LineWidth = 3;
%         end
%     else
%         warning('Speed 1 seems wrong, no 30 m/a contour')
%     end

    if(n == 1 || n < 0)
        ylabel('Northing [km]')
    end
    xlabel('Easting [km]')
    if(n == 4 || n < 0)
        c = colorbar;
        c.Label.String = 'Strain [yr^{-1}]';
    end
    title("\Delta" + "z = " + data1.thin_m + "m",'FontSize',ftsize)
%     caxis([1e-3 1e-1])     
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    f.ColorScale = 'linear';
    c.FontSize = ftsize;
    view(2)
    axis equal
    CT2 = cbrewer('seq','PuRd',256);
    colormap(ax, CT2);
%     xlim([-4.5e5 -2.5e5])
%     ylim([-5.75e5 -3.5e5])
end

