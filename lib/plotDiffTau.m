function[] = plotDiffTau(data1,data2,n,ax)
    ftsize = 20;
%     load Dawn.mat

    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        (data1.N_ef.*data1.tau_c(data1.xy(:,1), data1.xy(:,2),data1.u,data1.v)./norms([data1.u,data1.v],2,2))/1e3...
        - (data2.N_ef.*data2.tau_c(data2.xy(:,1), data2.xy(:,2),data2.u,data2.v)./norms([data2.u,data2.v],2,2))/1e3,...
        'edgecolor','none','facecolor','interp')
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
    CT = cmocean('tarn');
    colormap(ax, CT)
    caxis([-20 20])
    % axis off
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    c.FontSize = ftsize;
%     xlim([-1.56e6 -1.38e6])
%     ylim([-5.75e5 -3e5])
end
