clear
close all
load NoLakes_dhdt.mat
load grids/gridFlowRiseA02

[xx,yy] = meshgrid(xvec,yvec);

uu = measures_interp('speed',xx,yy);
uu = fillmissing(uu,'linear');
uu(isnan(dhdt)) = 0;

IN = inpolygon(xx,yy,pv(:,1),pv(:,2));
IN(uu < 30) = 0;

figure
scatter(xx(IN),yy(IN),[],dhdt(IN))
colorbar

figure
histogram(dhdt(IN))

meanDhDt = mean(mean(dhdt(IN)));
disp("Mean thinning is " + meanDhDt + " m/yr")
