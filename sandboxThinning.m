clear
rho_w = 1000;
rho   = 917;

%%
load ATL15_dhdt.mat
[Xi,Yi] = meshgrid(xvec,yvec);
spd = measures_interp('speed',Xi,Yi);

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


n_cases = 20;
randX = floor(rand(1,n_cases)*501);
randY = floor(rand(1,n_cases)*501);
figure(2)
clf
for i = 1:n_cases
    plot(squeeze(dhdt_pos(randX(i),randY(i),:)))
    hold on
    plot(squeeze(dhdt_neg(randX(i),randY(i),:)))
end

%% 
load Golledge21_GRL_T1_thick_22mar23_v2_Paul.mat

[Xi,Yi] = meshgrid(is2xvec,is2yvec);
spd = measures_interp('speed',Xi,Yi);

%%
figure(2)
clf
n = 6;
dn = floor(linspace(35,41,n));
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
dn = floor(linspace(1,41,n));
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
dn = floor(linspace(1,41,n));
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
HAF_end = (surf_interp(:,:,41) - bed_interp(:,:,41)) +  rho_w/rho*bed_interp(:,:,41);
surf(is2xvec,is2yvec,zeros(size(HAF_end)),HAF_end,'edgecolor', 'none')
hold on 
contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title("HAF at end")
view(2)
caxis([0 500])
colorbar
colormap parula

nexttile(2)
HAF_1 = (surf_interp(:,:,1) - bed_interp(:,:,1)) +  rho_w/rho*bed_interp(:,:,1);
surf(is2xvec,is2yvec,zeros(size(HAF_1)),HAF_1,'edgecolor', 'none')
hold on 
contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title("HAF at start")
view(2)
caxis([0 500])
colorbar
colormap parula

nexttile(3)
[Xi,Yi] = ndgrid(is2xvec,is2yvec);
HAF_bm = (bedmachine_interp('thickness',Xi,Yi)) +  rho_w/rho*bedmachine_interp('bed',Xi,Yi);
surf(is2xvec,is2yvec,zeros(size(HAF_bm)),HAF_bm,'edgecolor', 'none')
hold on 
contour(is2xvec,is2yvec,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(is2xvec,is2yvec,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title("HAF bedmachine")
view(2)
caxis([0 500])
colorbar
colormap parula