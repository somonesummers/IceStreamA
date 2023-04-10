clear
load ATL15_dhdt.mat

[Xi,Yi] = meshgrid(xvec,yvec);
spd = measures_interp('speed',Xi,Yi);

%%
figure(1)
clf
surf(xvec,yvec,zeros(size(dhdt)),dhdt,'edgecolor', 'none')
hold on 
contour(xvec,yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(xvec,yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(xvec,yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xvec,yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
view(2)
caxis([-1.5,1.5])
colorbar
colormap redblue

%% 
load Golledge21_GRL_T2_thick_22mar23_v2_Paul.mat

[Xi,Yi] = meshgrid(is2xvec,is2yvec);
spd = measures_interp('speed',Xi,Yi);

%%
figure(2)
clf
n = 6;
dn = floor(linspace(1,240,n));
tiledlayout(1,n+1,"TileSpacing","compact") 
for i = 1:n
    nexttile(i)
    contourf(is2xvec,is2yvec,surf_interp(:,:,dn(i)),'edgecolor', 'none')
    hold on 
    contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
    colorbar
    colormap parula
    title("Surface " + dn(i))
    caxis([0 3500])
    hold off    
end
    nexttile(n+1)
    contourf(is2xvec,is2yvec,bedmachine_interp('surface',Xi,Yi),'edgecolor', 'none')
    hold on 
    contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
    colorbar
    colormap parula
    title('Bedmachine')
    caxis([0 3500])
    hold off    

figure(3)
clf
n = 6;
dn = floor(linspace(1,240,n));
tiledlayout(1,n+1,"TileSpacing","compact") 
for i = 1:n
    nexttile(i)
    contourf(is2xvec,is2yvec,bed_interp(:,:,dn(i)),'edgecolor', 'none')
    hold on 
    contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
    colorbar
    colormap parula
    title("Bed " + dn(i))
    caxis([-1500 3000])
    hold off    
end
    nexttile(n+1)
    contourf(is2xvec,is2yvec,bedmachine_interp('bed',Xi,Yi),'edgecolor', 'none')
    hold on 
    contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
    colorbar
    colormap parula
    title('Bedmachine')
    caxis([-1500 3000])
    hold off    

 figure(4)
clf
n = 6;
dn = floor(linspace(1,240,n));
tiledlayout(1,n+1,"TileSpacing","compact") 
for i = 1:n
    nexttile(i)
    contourf(is2xvec,is2yvec,thick_interp(:,:,dn(i)),'edgecolor', 'none')
    hold on 
    contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
    colorbar
    colormap parula
    title("Thick " + dn(i))
    caxis([0 4000])
    hold off    
end
    nexttile(n+1)
    contourf(is2xvec,is2yvec,bedmachine_interp('thickness',Xi,Yi),'edgecolor', 'none')
    hold on 
    contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
    colorbar
    colormap parula
    title('Bedmachine')
    caxis([0 4000])
    hold off       
    
%%
figure(5)
clf
tiledlayout(2,2,"TileSpacing","compact") 

nexttile(1)
delta_bed = (bed_interp(:,:,end)-bed_interp(:,:,1));
surf(is2xvec,is2yvec,zeros(size(delta_bed)),delta_bed,'edgecolor', 'none')
hold on 
contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title("Change in bed height")
view(2)
% caxis([-1.5,1.5])
colorbar
colormap redblue

nexttile(2)
delta_surf = surf_interp(:,:,end)-surf_interp(:,:,1);
surf(is2xvec,is2yvec,zeros(size(delta_surf)),delta_surf,'edgecolor', 'none')
hold on 
contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title("Change in surface height")
view(2)
% caxis([-1.5,1.5])
colorbar
colormap redblue

nexttile(3)
BM_delta = thick(:,:,end)-bedmachine_interp('surface',is2xvec,is2xvec);
surf(is2xvec,is2yvec,zeros(size(BM_delta)),BM_delta,'edgecolor', 'none')
hold on 
contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title("Change in surface height end compared to today")
view(2)
% caxis([-1.5,1.5])
colorbar
colormap redblue