function [] = plotStress(data1,n,ax)
    ftsize = 20;
    tau_threshold = 100e3;
%     A = 3.5e-25;  %flow parameter pre-factor [Pa s^1/3] @-10C from cuffey 
%     A = 2.1e-2;  %flow parameter pre-factor [Pa s^1/3] @-15C from cuffey 
    A = 1.2e-25; % flow parameter pre-factor [Pa s^1/3] @ -20C from cuffey (~surf temp)
    nn = 3;
    tau = A^(-1/nn)*data1.ep_dot(data1.xy(:,1),data1.xy(:,2)).^(1/nn);
    
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
            zeros(size(data1.xy(:,1))),tau,'edgecolor','none','facecolor','interp');
    hold on
    
   
    if(max(tau) > tau_threshold)
        [~, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
            tau,[tau_threshold tau_threshold]);
        for j = 1:numel(H) %there can be multiple contours, so plot them all
            H(j).EdgeColor = rgb('lime green');
            H(j).LineStyle = ':';
            H(j).LineWidth = 2;
        end
    else
        warning('strain seems wrong, no strain threshold contour')
    end
    
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
        c.Label.String = 'Stress [Pa]';
    end
    title("\Delta" + "z = " + data1.thin_m + "m",'FontSize',ftsize)
    caxis([1e4 5e5])     
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    f.ColorScale = 'log';
    c.FontSize = ftsize;
    view(2)
    axis equal
    CT2 = cbrewer('seq','PuBuGn',256);
    colormap(ax, CT2);
%     xlim([-4.5e5 -2.5e5])
%     ylim([-5.75e5 -3.5e5])
end

