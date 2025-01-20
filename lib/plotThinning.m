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



