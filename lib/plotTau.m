function[] = plotTau(data1,n,ax)
    ftsize = 20;
    load Dawn.mat
    trisurf(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,zeros(size(data1.xy(:,1))),data1.N_ef.*data1.tau_c(data1.xy(:,1),...
        data1.xy(:,2),data1.u,data1.v)./norms([data1.u,data1.v],2,2)/1e3,...
       'edgecolor','none','facecolor','interp')
%     hold on
%     trisurf(t,xy(:,1),xy(:,2),-10e2*ones(size(xy(:,1))),...
%            'edgecolor','black','facecolor','none')
    caxis([0 250]);
    colormap(ax, (Cmap/255.0))
%     title("Basal Strength",'FontSize',ftsize+2)
    if(n == 1 || n < 0)
        ylabel('Northing [m]')
    end
    xlabel('Easting [m]')
    if(n == 4 || n < 0)
        c = colorbar;
        c.Label.String = 'Basal Strength [kPa]';
    end
    view(2)
    axis equal
    hold on 
    
%     if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
%         [~, H2] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
%             sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
%         for j = 1:numel(H2)
%         H2(j).EdgeColor = rgb('black');
%         H2(j).LineStyle = '--';
%         H2(j).LineWidth = 2;
%         end
%     else
%         warning('Speed 2 seems wrong, no 30 m/a contour')
%     end
    % axis off
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    c.FontSize = ftsize;
%     xlim([-1.56e6 -1.38e6])
%     ylim([-5.75e5 -3e5])
end

