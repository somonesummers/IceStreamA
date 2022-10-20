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

xbox = val.features.geometry.coordinates(:,:,1);
ybox = val.features.geometry.coordinates(:,:,2);




icey = cbrewer('div','BrBG',48);
rho = 917;
rho_w = 1000;
g = 9.81;
B = 1.6e8; % A = 2.4e-25 Pa^(-3) s^(-1)
A = 2.4e-25;
overgrab = 0;
xmax =  -2.5e5;
xmin =  -5.0e5;
ymax =  -3.5e5;
ymin =  -6.0e5;

dx = 2e3;
smth = 10e3;
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
e_eff = sqrt(.5*(ux.^2 + vy.^2) + (.5*(uy + vx)).^2);
[e_effx, e_effy] = gradient(e_eff.^(1/3-1),dx,dx);

% Deviatoric Stress
prefactor = B * e_eff.^(1/3 - 1);

% Driving Stress
Tdx = -rho * g * h .* sx.^2;
Tdy = -rho * g * h .* sy.^2;
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
        dr(j,i) = -(vv(1)*sx(j,i) + vv(2)*sy(j,i))* rho * g * h(j,i); %Driving Force
        
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
% u_int = dr^3 * H * A;

u_int = 2 / 4 *abs(bed).^3 .* h * A * 3.154e7;
%%
figure(1)
clf
ax(1) = subplot(121);
p = surf(Xi,Yi,zeros(size(ss)),log10(spd2));
hold on 
utmp = cos(angs);
vtmp = sin(angs);
% quiver(Xi,Yi,utmp,vtmp)
title('Ice Speed')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Log_{10} Speed [m/yr]';
ax(2) = subplot(122);
% p = surf(Xi,Yi,zeros(size(ss)),b_raw);
p = surf(Xi,Yi,zeros(size(b_raw)),b_raw);
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Bed elevation')
set(p, 'edgecolor', 'none');
colormap(ax(2),icey);
% caxis([-2000 500])
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Bed Elevation [m]';

%% Hydropotential
figure(2)
clf
tiledlayout(1,2, 'Padding', 'tight', 'TileSpacing', 'tight');
ax1 = nexttile(1);
p = surf(Xi,Yi,zeros(size(phi_raw)),phi_raw);
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)

contour(xi,yi,phi_raw, [0:10:700] ,'color',rgb('purple'),'LineWidth',1,'HandleVisibility','off')
contour(xi,yi,phi_raw, [0:50:700] ,'color',rgb('purple'),'LineWidth',2,'HandleVisibility','off')
title('Hydro Potential (hydroequalib assumption)')
set(p, 'edgecolor', 'none');
colormap(ax1,cbrewer('seq','BuPu',256));
caxis([0 700])
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = '\Phi [?]';

ax2 = nexttile(2);
[phi_x, phi_y] = gradient(phi,dx,dx);
phi_x_n = phi_x./sqrt(phi_x.^2 + phi_y.^2);
phi_y_n = phi_y./sqrt(phi_x.^2 + phi_y.^2);
[b_x, b_y] = gradient(b,dx,dx);
sp = 3;
p = surf(Xi,Yi,zeros(size(phi_raw)),sqrt(phi_x.^2 + phi_y.^2));
hold on
set(p, 'edgecolor', 'none');
quiver(xi(1:sp:end),yi(1:sp:end),-phi_x_n(1:sp:end,1:sp:end),-phi_y_n(1:sp:end,1:sp:end));
% quiver(xi(1:sp:end),yi(1:sp:end),-b_x(1:sp:end,1:sp:end),-b_y(1:sp:end,1:sp:end));
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
hold off
colorbar
colormap(ax2,flipud(pink))
caxis([0 .05])
title('Gradients')
legend('|\nabla \Phi|','\nabla \Phi')
view(2)
axis equal
setFontSize(16);

%% 
figure(3)
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
% figure(3)
% clf
% sgtitle('Force Budget Fraction (Positive is Along Flow)')
% colormap redblue
% caxis([-1e5 1e5])
% subplot(221)
% p = surf(Xi,Yi,zeros(size(ss)),dr);
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Driving Force')
% allfig(p)
% 
% subplot(222)
% p = surf(Xi,Yi,zeros(size(ss)),lon./dr);
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Longitudinal Stresses')
% allfig(p)
% caxis([-1, 1])
% c = colorbar;
% c.Label.String = '[%]';
% 
% subplot(223)
% p = surf(Xi,Yi,zeros(size(ss)),lat./dr);
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Lateral Stresses')
% allfig(p)
% caxis([-1, 1])
% c = colorbar;
% c.Label.String = '[%]';
% 
% subplot(224)
% p = surf(Xi,Yi,zeros(size(ss)),bed./dr);
% title('Bed Drag')
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% allfig(p)
% caxis([-1, 1])
% c = colorbar;
% c.Label.String = '[%]';
% 
% setFontSize(16)
%% 

%%
figure(5)
clf
colormap redblue
p = surf(Xi,Yi,zeros(size(ss)),lat);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Lateral Stresses')
bedmachine('gl','k','linewidth',2)
allfig2(p,lat)

%%

figure(6)
clf
ax1 = axes;
r = 4;
surf(ax1,Xi(1:r:end,1:r:end),Yi(1:r:end,1:r:end),b_raw(1:r:end,1:r:end),...
    b_raw(1:r:end,1:r:end),'facealpha',[.9],'edgecolor', 'k');
colormap(copper)
% lighting gouraud
c = colorbar('south');
c.Label.String = 'Bed Elevation [m ASL]';
hold on
ax2 = axes;
p = surf(ax2,Xi,Yi,sf_raw,(spd2),'facealpha',0.9);
title('Elevation with Bed')
set(p, 'edgecolor', 'none');
colormap(parula);
% caxis([-2000 2000])
% axis equal
setFontSize(16);
c = colorbar('north');
c.Label.String = 'Ice Surface Speed [m/yr]';
f = gca;
f.ColorScale = 'log';
set(p, 'edgecolor', 'none');

hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'copper')
colormap(ax2,'parula')
%%

figure(7)
clf

p = surf(Xi,Yi,sf_raw,(spd2),'facealpha',1);
hold on
contour3(Xi,Yi,sf_raw,[-1000:100:2000],'-','color',rgb('gray'))
contour3(Xi,Yi,sf_raw,[-1000:500:2000],'k-')
title('Elevation')
set(p, 'edgecolor', 'none');
colormap(parula);
% caxis([-2000 2000])
% axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Ice Surface Speed [m/yr]';
f = gca;
f.ColorScale = 'log';
set(p, 'edgecolor', 'none');

%%
figure(8)
clf
subplot(121)
p = surf(Xi,Yi,zeros(size(ss)),u_int);
title('Vertical Creep (Idealized)')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,u_int)
caxis([0 100])
c = colorbar;
c.Label.String = '[m/yr]';


subplot(122)
p = surf(Xi,Yi,zeros(size(ss)),(abs(u_int)./spd2*100));
title('Internal Creep Factor')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,u_int./spd2*100)
c = colorbar;
caxis([0,100]);
c.Label.String = '[%]';
f = gca;
f.ColorScale = 'linear';
%%
figure(9)
clf
p = surf(Xi,Yi,zeros(size(ss)),h);
title('Thickness')
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,h)
c = colorbar;
c.Label.String = '[m]';
caxis([0 3000])
%%
figure(10)
clf
p = contourf(Xi,Yi,imgaussfilt(b_raw,3),10);
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Bed elevation')
colormap(ax(2),icey);
caxis([-2000 500])
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Bed Elevation [m]';
%%
figure
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
quiver(xi(1:8:end),yi(1:8:end),-sx_hat(1:8:end,1:8:end),-sy_hat(1:8:end,1:8:end),.5,'k')
title('Surface Slope')
bedmachine('gl','k','linewidth',2)
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
set(p, 'edgecolor', 'none');
caxis([1e-6 1e-2])
%%
figure
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

% figure
% clf
% sgtitle('Input sources isnan')
% subplot(221)
% surf(Xi,Yi,zeros(size(ss)),int16(isnan(b)),'edgecolor','none');
% view(2)
% colorbar
% title('bed')
% 
% subplot(222)
% surf(Xi,Yi,zeros(size(ss)),int16(isnan(sf)),'edgecolor','none');
% view(2)
% colorbar
% title('surface')
% 
% 
% subplot(223)
% surf(Xi,Yi,zeros(size(ss)),int16(isnan(u)),'edgecolor','none');
% view(2)
% colorbar
% title('u')
% 
% subplot(224)
% surf(Xi,Yi,zeros(size(ss)),int16(isnan(v)),'edgecolor','none');
% view(2)
% colorbar
% title('v')

figure
clf

colormap redblue
p = surf(Xi,Yi,zeros(size(ss)),abs(dr)-abs(bed),'edgecolor','none');
view(2)
hold on
colorbar
caxis([-1.5 1.5]*1e5)
title('abs(Driving Force) - abs(Drag)')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','c-','linewidth',2)



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

