clear

fname = ('ice-stream-a-domain.geojson');
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

xbox = val.features.geometry.coordinates(:,:,1);

ybox = val.features.geometry.coordinates(:,:,2);

[bed,R1] = readgeoraster('bedmachineAntarctica_v02_bed_ross.tif');
[surface,R1] = readgeoraster('bedmachineAntarctica_v02_surface_ross.tif');
[xvel,R2] = readgeoraster('antarctic_ice_vel_phase_v01_vx_ross.tif');
[yvel,R2] = readgeoraster('antarctic_ice_vel_phase_v01_vy_ross.tif');

% R1 = flipud(R1);
% R2 = flipud(R2);
bed = flipud(bed);
surface = flipud(surface);
xvel = flipud(xvel);
yvel = flipud(yvel);

xx1 = R1.XWorldLimits(1):R1.CellExtentInWorldX:R1.XWorldLimits(2);
yy1 = R1.YWorldLimits(1):R1.CellExtentInWorldY:R1.YWorldLimits(2);
[XX1, YY1] = meshgrid(xx1(2:end),yy1(2:end));

xx2 = R2.XWorldLimits(1):R2.CellExtentInWorldX:R2.XWorldLimits(2);
yy2 = R2.YWorldLimits(1):R2.CellExtentInWorldY:R2.YWorldLimits(2);
[XX2, YY2] = meshgrid(xx2(2:end),yy2(2:end));

%%
% figure(1)
% clf
% surf(XX1,YY1,bed,'edgecolor','none')
% hold on
% plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
% view(2)
% title('Bed Height [m ASL]')
% colorbar
% 
% figure(2)
% clf
% surf(XX2,YY2,sqrt(xvel.^2 + yvel.^2),'edgecolor','none')
% hold on
% plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
% title('Surface Speed [m/yr]')
% view(2)
% colorbar

b_err = bedmachine_interp('errbed',XX1,YY1);
method = bedmachine_interp('source',XX1,YY1);

figure
clf
surf(XX2,YY2,zeros(size(xvel)),sqrt(xvel.^2 + yvel.^2),'edgecolor','none')
hold on
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[30,100,300],'k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[10,10],':k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[1000,1000],'k','linewidth',2)
view(2)
title('Speed with Speed Contours')
xlim([-6.5, 0]*1e5);
ylim([-7.2, -1.5]*1e5);
set(gca,'ColorScale','log')
caxis([1 2000])
colorbar


figure
clf
surf(XX1,YY1,zeros(size(bed)),surface - bed,'edgecolor','none')
hold on
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[30,100,300],'k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[10,10],':k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[1000,1000],'k','linewidth',2)
view(2)
title('Thickness [m] with Speed Contours, BedMachine')
xlim([-6.5, 0]*1e5);
ylim([-7.2, -1.5]*1e5);
caxis([0 4000])
colorbar
%%
figure
clf
surf(XX1,YY1,zeros(size(bed)),b_err,'edgecolor','none')
hold on
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[30,100,300],'k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[10,10],':k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[1000,1000],'k','linewidth',2)
view(2)
title('Bed Error [m] with Speed Contours, BedMachine')
xlim([-6.5, 0]*1e5);
ylim([-7.2, -1.5]*1e5);
colorbar
caxis([200 600])

%%
figure
clf
surf(XX1,YY1,zeros(size(bed)),method,'edgecolor','none')
hold on
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[30,100,300],'k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[10,10],':k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[1000,1000],'k','linewidth',2)
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
view(2)
title('Method with Speed Contours, BedMachine')
%'1 = REMA/IBCSO, 2 = Mass conservation, 3 = interpolation, 4 = hydrostatic, 5 = Kriging/streamline?, 6 = gravity inversion');
xlim([-6.5, 0]*1e5);
ylim([-7.2, -1.5]*1e5);
colorbar
colormap jet
caxis([0 6])

%%
bedDiff = bedmap2_interp(XX1,YY1,'bed')-bed;

figure
clf
surf(XX1,YY1,zeros(size(bed)),bedDiff,'edgecolor','none')
hold on
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[30,100,300],'k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[10,10],':k')
contour3(XX2,YY2,sqrt(xvel.^2 + yvel.^2),[1000,1000],'k','linewidth',2)
plot3(xbox,ybox,2000*ones(size(xbox)),'k','linewidth',3)
view(2)
title('BedMap - BedMachine Bed Elevation [m]')
xlim([-6.5, 0]*1e5);
ylim([-7.2, -1.5]*1e5);
colorbar
colormap redblue
caxis([-1500 1500])