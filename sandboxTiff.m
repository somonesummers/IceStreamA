clear
% data = Tiff('IS2_smoothed_no_lakes.tif','r');
% values = imread(imread(data));
values = double(flipud(imread('IS2_smoothed_no_lakes.tif')));


load('ATL15_dhdt.mat','xvec','yvec','dhdt')
PS = load('PS_dhdt.mat');

values(isnan(dhdt)) = nan;

% mean(mean(values,'omitnan'))
% mean(mean(dhdt,'omitnan'))
% 
% scale = mean(mean(dhdt,'omitnan'))/mean(mean(values,'omitnan'));

values = values/2.7491; %convert to m/yr

figure(1)
clf
surf(xvec/1e3,yvec/1e3,values,'edgecolor','none')
view(2)
colormap(redblue)
caxis([-1.5 1.5])
c = colorbar;
c.Label.String = 'Height Change [m/yr]';
xlabel("Easting [km]")
ylabel("Northing [km]")
title('No Lakes')
setFontSize(18)

figure(2)
clf
surf(xvec/1e3,yvec/1e3,dhdt,'edgecolor','none')
view(2)
colormap(redblue)
caxis([-1.5 1.5])
c = colorbar;
c.Label.String = 'Height Change [m/yr]';
xlabel("Easting [km]")
ylabel("Northing [km]")
title('ALT15 dhdt')
setFontSize(18)


figure(3)
clf
surf(xvec/1e3,yvec/1e3,dhdt-values,'edgecolor','none')
view(2)
colormap(redblue)
caxis([-1.5 1.5])
c = colorbar;
c.Label.String = 'Difference [m/yr]';
xlabel("Easting [km]")
ylabel("Northing [km]")
title('ALT - NoLakes')
setFontSize(18)

dhdt = values;
save("NoLakes_dhdt",'xvec','yvec','dhdt');
