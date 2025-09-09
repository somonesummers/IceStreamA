% Mapping strain and driving force for ice stream A methods of Van Der Veen 1989
clear; close all;
addpath lib
load Dawn.mat

fname = ('ice-stream-a/ice-stream-a-domain.geojson');
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

xbox = [-4.9500  -3.6000   -2.7000   -4.1000   -4.9500]*1e5;
ybox = [-3.8000  -3.5300   -5.2250   -5.6700   -3.8000]*1e5;

xmax =  max(xbox);
xmin = min(xbox);
ymax =  max(ybox);
ymin =  min(ybox);

pv = [xbox; ybox]';




icey = cbrewer('div','BrBG',48);
rho = 917;
rho_w = 1000;
g = 9.81;
% B = 1.419e8; % A = 3.5e-25 Pa^(-3) s^(-1)
A = 3.5e-25;
B = A^(-1/3);
overgrab = 50;
% xmax =  -2.5e5;
% xmin =  -5.0e5;
% ymax =  -3.5e5;
% ymin =  -6.0e5;

dx = 1e3;
smth = 4e3;
xi = xmin-dx*overgrab:dx:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% Masking
% msk = scatteredInterpolant(xy(:,1),xy(:,2),ones(size(xy(:,1))),'nearest','none');
% mask = msk(Xi,Yi);
% mask(isnan(mask)) = 0;
 mask = ones(size(Xi)); %make mask all ones for now

%Raw fields
b_raw =  bedmachine_interp('bed',Xi,Yi);
sf_raw =  bedmachine_interp('surface',Xi,Yi);
% b_raw =  bedmap2_interp(Xi,Yi,'bed');
% sf_raw =  bedmap2_interp(Xi,Yi,'surface');
phi_raw = rho/rho_w.*sf_raw + (rho_w-rho)/rho_w.*b_raw;

[u, v] = measures_interp('velocity',Xi,Yi);

%% Cleaning NANs in velocity [do cautiously]
valid_u     = ~isnan(u);
u_interp = scatteredInterpolant(Xi(valid_u),Yi(valid_u),u(valid_u),'natural');
u(~valid_u) = u_interp(Xi(~valid_u),Yi(~valid_u));
valid_v     = ~isnan(v);
v_interp = scatteredInterpolant(Xi(valid_v),Yi(valid_v),v(valid_v),'natural');
v(~valid_v) = v_interp(Xi(~valid_v),Yi(~valid_v));


%% Smoothing/Masking
% Gauss
b = imgaussfilt(b_raw,2);
u = imgaussfilt(u,smth/dx) / 3.154e7;
v = imgaussfilt(v,smth/dx) / 3.154e7;
sf = imgaussfilt(sf_raw,smth/dx);
phi = imgaussfilt(phi_raw,2);
% None
% b = b_raw;
% u = u/ 3.154e7;
% v = v / 3.154e7;
% sf = sf_raw;
% phi = phi_raw;

spd = sqrt(u.^2 + v.^2);
% mean filtering
% hh = ones(floor(smth/dx))/(floor(smth/dx).^2);
% u = filter2(hh,u) / 3.154e7;
% v = filter2(hh,v) / 3.154e7;
% surf = filter2(hh,surf);

sf(sf < b) = b(sf < b);
spd2  = measures_interp('speed',Xi,Yi);

%Gradients
h = sf-b;
[ux ,  uy] = gradient(u,dx,dx);
[vx ,  vy] = gradient(v,dx,dx);
[sx ,  sy] = gradient(sf,dx,dx);
[spdx ,  spdy] = gradient(spd,dx,dx);
[spdxx, spdxy] = gradient(spdx,dx,dx);
[spdyx, spdyy] = gradient(spdy,dx,dx);


% Effective Strain
e_eff = sqrt(.5*(ux.^2 + vy.^2) + (.5*(ux + uy).^2) + (.25*(uy + vx)).^2);
[e_effx, e_effy] = gradient(e_eff.^(1/3-1),dx,dx);

% Deviatoric Stress
prefactor = B * e_eff.^(1/3 - 1);

% Driving Stress
Tdx = -rho * g * h .* sx;
Tdy = -rho * g * h .* sy;
Td  = sqrt(Tdx.^2 +  Tdy.^2);

ss   = zeros(size(u));
angs = zeros(size(u));
dr   = zeros(size(u));
lon  = zeros(size(u));
lat  = zeros(size(u));
bed  = zeros(size(u));
h = sf-b;
for i = 2:length(xi)-1
    for j = 2:length(yi)-1
        ui = u(j,i);
        vi = v(j,i);
        ang = atan(vi/ui);
        if(ui < 0) 
            ang = ang + pi;
        end
        vv = [cos(ang), sin(ang)]; %Direction Vectors along flow
        vv_t = [-sin(ang), cos(ang)];%Direction Vectors Perp to flow
        R = [ cos(ang) sin(-ang) ; sin(ang) cos(ang) ];
        ss(j,i) = max([ui , vi] * R); %x speed in new ref frame
        % dr(j,i) = -(vv(1)*sx(j,i) + vv(2)*sy(j,i))* rho * g * h(j,i); %Driving Force
        dr(j,i) = Td(j,i);
        lon(j,i) =  2*B*((vv(1)*sx(j,i) + vv(2)*sy(j,i)) .* e_eff(j,i).^(1/3-1) .* (vv(1)*spdx(j,i) + vv(2)*spdy(j,i))...
                    + h(j,i) .* (vv(1)*e_effx(j,i) + vv(2)*e_effy(j,i)) .* (vv(1)*spdx(j,i) + vv(2)*spdy(j,i))...
                    + h(j,i) .* e_eff(j,i).^(1/3-1) .* (spdxx(j,i).*vv(1).^2 + spdxy(j,i).*vv(1).*vv(2) + spdyx(j,i).*vv(1).*vv(2) + spdyy(j,i).*vv(2).^2));
                
        lat(j,i) =  2*B*((vv_t(1)*sx(j,i) + vv_t(2)*sy(j,i)) .* e_eff(j,i).^(1/3-1) .* (vv_t(1)*spdx(j,i) + vv_t(2)*spdy(j,i))...
                    + h(j,i) .* (vv_t(1)*e_effx(j,i) + vv_t(2)*e_effy(j,i)) .* (vv_t(1)*spdx(j,i) + vv_t(2)*spdy(j,i))...
                    + h(j,i) .* e_eff(j,i).^(1/3-1) .* (spdxx(j,i).*vv_t(1).^2 + spdxy(j,i).*vv_t(1).*vv_t(2) + spdyx(j,i).*vv_t(1).*vv_t(2) + spdyy(j,i).*vv_t(2).^2));
        angs(j,i) = ang;
        
        bed(j,i) = dr(j,i) + lat(j,i) + lon(j,i);
    end
end

% Masking
dr = dr .* mask;
lat = lat .* mask;
lon = lon .* mask;
bed = bed .* mask;
spd2 = spd2 .* mask;

% Internal Deformation expected over a locked bed [m/yr]
% u_int = 2*A / (n+1) dr^3 *H ; %Cuffey eq 8.35

u_int = 2 / 4 *abs(dr).^3 .* h * A * 3.154e7;
max_vert_shear = A * dr.^3 * 3.154e7; %[1/year]
%log stain rates using face formula
[elon, etrans, eshear, x_strain, y_strain] = strainRates([Xi(:),Yi(:)],u(:),v(:),h(:)); 
%more basic approach
alpha = atan(v./u);
e_xy = .5 * (uy + vx);
e_shr = (vy-ux).*cos(alpha).*sin(alpha) + e_xy.*(cos(alpha).^2 - (sin(alpha).^2));
%
%%
figure(6)
tiledlayout(1,3, 'Padding', 'tight', 'TileSpacing', 'tight');
clf
sgtitle('Strain Estimates')
colormap parula
ax = nexttile(1);
p = surf(Xi,Yi,zeros(size(ss)),max_vert_shear);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
view(2)
clim([2e-5 2e-3]);
set(p, 'edgecolor', 'none');
xlabel('Easting')
ylabel('Northing')
title('Estimated Vertical Shear Strain Maximum')
axis equal
c = colorbar;
c.Label.String = '[1/year ]';

ax = nexttile(2);
% p = surf(x_strain,y_strain,zeros(size(eshear)),abs(eshear)*3.154e7);
p = surf(xi,yi,zeros(size(e_shr)),abs(e_shr)*3.154e7);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
view(2)
clim([2e-5 2e-3]);
set(p, 'edgecolor', 'none');
xlabel('Easting')
ylabel('Northing')
title('Estimated Lateral Shear Strain from Obs')
axis equal
c = colorbar;
c.Label.String = '[1/year ]';

ax = nexttile(3);
p = surf(xi,yi,zeros(size(e_shr)),abs(e_shr)*3.154e7./max_vert_shear);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
view(2)
clim([0 1]);
set(p, 'edgecolor', 'none');
colormap(ax,redblue)
xlabel('Easting')
ylabel('Northing')
title('Ratio of Horizontal to Vertical Shear')
axis equal
c = colorbar;
c.Label.String = '[ ]';

%% 
figure(1)
clf
sgtitle('Force Budget (Positive is Along Flow)')
colormap redblue
caxis([-1e5 1e5])
subplot(221)
p = surf(Xi,Yi,zeros(size(ss)),dr);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
title('Driving Force')
allfig2(p,dr)

subplot(222)
p = surf(Xi,Yi,zeros(size(ss)),lon);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
title('Longitudinal Stresses')
allfig2(p,lon)

subplot(223)
p = surf(Xi,Yi,zeros(size(ss)),lat);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
title('Lateral Stresses')
allfig2(p,lat)

subplot(224)
p = surf(Xi,Yi,zeros(size(ss)),bed);
title('Bed Drag')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
contour(xi,yi,h, [100, 100] , 'g-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
allfig2(p,bed)

% bed(h<100) = 200e3;
% save("bedDragDx" + string(dx) + "smth" + string(smth) +"Capped.mat",'Xi','Yi','bed')

save("bedDragDx" + string(dx) + "smth" + string(smth) +".mat",'Xi','Yi','bed')

setFontSize(16)

%%
figure(2)
clf
spd2 = fillmissing(spd2,'linear');
p = surf(Xi/1e3,Yi/1e3,zeros(size(ss)),(spd2 - abs(u_int)),'edgecolor','none');
% title('Estimated Basal Sliding')
view(2)
hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi/1e3,yi/1e3,spd2, [5,15,25] , 'k-','HandleVisibility','off')
contour(xi/1e3,yi/1e3,spd2, [10,20,30] , 'k-','LineWidth',2)
% allfig2(p,u_int./spd2*100)
c = colorbar;
axLocal = gca;
colormap(cbrewer('seq','Greens',64));
caxis([0,30]);
c.Label.String = 'Est Basal Speed [m/yr]';
f = gca;
f.ColorScale = 'linear';
axis equal
setFontSize(24)
% axis off
xlabel("Easting [km]")
ylabel("Northing [km]")
xlim([-4.7347   -2.8804]*1e2);
ylim([-5.5556   -3.7640]*1e2);
savePng('figs/basalSliding')

%%
figure(3)
p = surf(Xi,Yi,zeros(size(ss)),sqrt(sx.^2+sy.^2));
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
[C,hh] = contour(xi,yi,sqrt(sx.^2+sy.^2), [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1], 'r-','HandleVisibility','off');  
clabel(C,hh)
hold on
sx_hat = sx./sqrt(sx.^2+sy.^2);
sy_hat = sy./sqrt(sx.^2+sy.^2);
u_hat = u./sqrt(u.^2 + v.^2);
v_hat = v./sqrt(u.^2 + v.^2);
quiver(xi(1:8:end),yi(1:8:end),-sx_hat(1:8:end,1:8:end),-sy_hat(1:8:end,1:8:end),.5,'k')
quiver(xi(1:8:end),yi(1:8:end),u_hat(1:8:end,1:8:end),v_hat(1:8:end,1:8:end),.5,'r')
title('Surface Slope (k = surface slope, r = velocity)')
bedmachine('gl','k','linewidth',2)
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
set(p, 'edgecolor', 'none');
caxis([1e-6 1e-2])
%%
figure(4)
p = surf(Xi,Yi,zeros(size(ss)),measures_interp('err',Xi,Yi)./measures_interp('speed',Xi,Yi),...
    'edgecolor', 'none');
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
view(2)
f = gca;
caxis([0 1])
colorbar
title('Percent error in speed')

figure(5)
clf
sgtitle('Fractional Contributions Lat/Driving (positive along flow)')
colormap redblue
p = surf(Xi,Yi,zeros(size(ss)),abs(lat)./abs(dr));
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)
view(2)
caxis([0 2]);
set(p, 'edgecolor', 'none');
xlabel('Easting')
ylabel('Northing')
axis equal
c = colorbar;
c.Label.String = '[ ]';


function [] = allfig(p)
set(p, 'edgecolor', 'none');
view(2)
c = colorbar;
c.Label.String = '[Pa]';
% caxis([-1.5e5, 1.5e5])
xlabel('Easting')
ylabel('Northing')
axis equal
end

function [] = allfig2(p,z)
set(p, 'edgecolor', 'none');
view(2)
c = colorbar;
c.Label.String = '[Pa]';
si = std(abs(z),0,'all','omitnan');
me = mean(mean(abs(z),'omitnan'));
caxis([-(me+3*si), (me+3*si)])
% caxis([-1.5e5, 1.5e5])
xlabel('Easting')
ylabel('Northing')
axis equal
end

