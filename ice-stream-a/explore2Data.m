clear;
addpath ../lib

fname = ('ice-stream-a-domain.geojson');
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

xbox = val.features.geometry.coordinates(:,:,1);
ybox = val.features.geometry.coordinates(:,:,2);

icey = cbrewer('div','BrBG',48);
rho = 917;
rho_w = 1000;
g = 9.81;
B = 1.6e8; % A = 2.4e-25 Pa^(-3) s^(-1)
A = 2.4e-25;
overgrab = 20;
xmax =  max(xbox);
xmin = min(xbox);
ymax =  max(ybox);
ymin =  min(ybox);

dx = 1e3;
smth = 5e3;
xi = xmin-dx*overgrab:dx:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% Masking
% msk = scatteredInterpolant(xy(:,1),xy(:,2),ones(size(xy(:,1))),'nearest','none');
% mask = msk(Xi,Yi);
% mask(isnan(mask)) = 0;
 mask = ones(size(Xi)); %make mask all ones for now

%Raw fields
% b_raw =  bedmachine_interp('bed',Xi,Yi);
% sf_raw =  bedmachine_interp('surface',Xi,Yi);
b_raw =  bedmap2_interp(Xi,Yi,'bed');
sf_raw =  bedmap2_interp(Xi,Yi,'surface');

[u, v] = measures_interp('velocity',Xi,Yi);

%% Smoothing/Masking
% Gauss
b = imgaussfilt(b_raw,2);
u = imgaussfilt(u,smth/dx) / 3.154e7;
v = imgaussfilt(v,smth/dx) / 3.154e7;
spd = sqrt(u.^2 + v.^2);
sf = imgaussfilt(sf_raw,smth/dx);


% mean filtering
% hh = ones(floor(smth/dx))/(floor(smth/dx).^2);
% u = filter2(hh,u) / 3.154e7;
% v = filter2(hh,v) / 3.154e7;
% surf = filter2(hh,surf);

sf(sf < b) = b(sf < b);
spd2  = measures_interp('speed',Xi,Yi);

%%
b_err =  bedmap2_interp(Xi,Yi,'beduncertainty');
figure
clf
surf(Xi,Yi,zeros(size(spd2)),sf_raw - b_raw,'edgecolor','none')
hold on
contour3(Xi,Yi,spd2,[30,100,300],'k')
contour3(Xi,Yi,spd2,[10,10],':k')
contour3(Xi,Yi,spd2,[1000,1000],'k','linewidth',2)
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
view(2)
title('Ice Thickness [m] with Speed Contours, BedMap2')
colorbar
caxis([0 4000])

figure
clf
surf(Xi,Yi,zeros(size(spd2)),b_err,'edgecolor','none')
hold on
contour3(Xi,Yi,spd2,[30,100,300],'k')
contour3(Xi,Yi,spd2,[10,10],':k')
contour3(Xi,Yi,spd2,[1000,1000],'k','linewidth',2)
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
view(2)
title('Bed Uncert [m] with Speed Contours, BedMap2')
colorbar
caxis([200 600])

