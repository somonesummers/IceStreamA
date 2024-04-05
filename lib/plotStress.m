function [] = plotStress(data1,n,ax)
    ftsize = 20;
    tau_threshold = 200e3;
    tau_threshold_med = 150e3;
    tau_threshold_low = 100e3;
%     A = 3.5e-25;  %flow parameter pre-factor [Pa s^1/3] @-10C from cuffey 
%     A = 2.1e-2;  %flow parameter pre-factor [Pa s^1/3] @-15C from cuffey 
    A = 1.2e-25; % flow parameter pre-factor [Pa s^1/3] @ -20C from cuffey (~surf temp)
    nn = 3;
    
    tau = sqrt(3)*A^(-1/nn)*data1.ep_dot(data1.xy(:,1),data1.xy(:,2)).^(1/nn); 
    %sqrt 3 for von mesis stress from Mercenier 2019, Vaughan 1993.
    
    trisurf(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            zeros(size(tau)),tau,'edgecolor','none','facecolor','interp');
    hold on
    
   
    if(max(tau) > tau_threshold)
        [~, H] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            tau,[tau_threshold tau_threshold]);
        for j = 1:numel(H) %there can be multiple contours, so plot them all
            H(j).EdgeColor = rgb('lime green');
            H(j).LineStyle = ':';
            H(j).LineWidth = 2;
        end
        [~, H] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            tau,[tau_threshold_med tau_threshold_med]);
        for j = 1:numel(H) %there can be multiple contours, so plot them all
            H(j).EdgeColor = rgb('robin egg blue');
            H(j).LineStyle = ':';
            H(j).LineWidth = 2;
        end
        [~, H] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
            tau,[tau_threshold_low tau_threshold_low]);
        for j = 1:numel(H) %there can be multiple contours, so plot them all
            H(j).EdgeColor = rgb('gray');
            H(j).LineStyle = ':';
            H(j).LineWidth = 2;
        end
        
    else
        warning('strain seems wrong, no strain threshold contour')
    end
    
    if(max(sqrt(data1.u.^2 + data1.v.^2)*3.154E7) > 30)
        [~, H] = tricontour(data1.t,data1.xy(:,1)/1e3,data1.xy(:,2)/1e3,...
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
    xlabel('Easting [km]')
    if(n == 4 || n < 0)
        c = colorbar;
        c.Label.String = 'Stress [Pa]';
    end
    title("\Delta" + "z = " + data1.thin_m + "m",'FontSize',ftsize)
    caxis([0 2*tau_threshold]);
    [cmin, cmax] = caxis();     
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    f.ColorScale = 'linear';
    c.FontSize = ftsize;
    view(2)
    axis equal
    %force abrubt color jump at Tau threshold on log scale
    CT1 = cbrewer('seq','PuBu',256);
    CT2 = cbrewer('seq','BuGn',256);
    clin = linspace(cmin,cmax,256); %swap clin for clog if you want linear
%     clog = logspace(log10(cmin),log10(cmax),256);
    [~,ind] = min(abs(clin-tau_threshold));
    CTnew = [CT1(1:ind,:); CT2(1+ind:end,:)];
%     [~,ind_m] = min(abs(clin-tau_threshold_med));
%     [~,ind_l] = min(abs(clin-tau_threshold_low));
%     CTnew(ind_m,:) = rgb('robin egg blue');
%     CTnew(ind_l,:) = rgb('gray');
    
    colormap(ax, CTnew);
end

