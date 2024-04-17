function [] = plotLogStrain(data1,n,ax)
% plotLogStrain(data1,n,ax) n = 1,2,3 : lat,trans,shear
    
    ftsize = 20;
    
    [elon, etrans, eshear, Xi, Yi] = strainRates(data1.xy,data1.u,data1.v,data1.h);
    inGrid = inpolygon(Xi,Yi,data1.pv(:,1),data1.pv(:,2));
    elon(~inGrid) = nan;
    etrans(~inGrid) = nan;
    eshear(~inGrid) = nan;
    if(n==1)
        surf(Xi/1e3,Yi/1e3,zeros(size(elon)),elon*3.154e7,'edgecolor','none')
    elseif(n==2)
        surf(Xi/1e3,Yi/1e3,zeros(size(elon)),etrans*3.154e7,'edgecolor','none')
    else
        surf(Xi/1e3,Yi/1e3,zeros(size(elon)),eshear*3.154e7,'edgecolor','none')
    end
    view(2)
    hold on
    
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
    caxis([-6e-2 6e-2])     
    view(2)
    f = gca;
    f.XAxis.FontSize = ftsize-2;
    f.YAxis.FontSize = ftsize-2;
    f.ColorScale = 'linear';
    c = colorbar;
    c.FontSize = ftsize;
    view(2)
    axis equal
    CT2 = cbrewer('div','PuOr',256);
    colormap(ax, CT2);





end