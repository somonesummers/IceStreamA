clear
load grids/gridFlowRiseC035.mat

overgrab = 20;
xxx = xmin-dx*overgrab:dx/2:xmax+dx*overgrab;
yyy = ymin-dx*overgrab:dx/2:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xxx,yyy);

% % % Raw fields
% bm2_b =  bedmap2_interp(Xi,Yi,'bed');
% bm2_s =  bedmap2_interp(Xi,Yi,'surface');
% 
% bm_b =  bedmachine_interp('bed',Xi,Yi);
% bm_s =  bedmachine_interp('surface',Xi,Yi);

[u,v] = measures_interp('velocity',Xi,Yi);
spd2 = measures_interp('speed',Xi,Yi);


%%
figure(1)
clf
quiver(Xi,Yi,u,v)
hold on
contour(xxx,yyy,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xxx,yyy,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xxx,yyy,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xxx,yyy,spd2, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
title('Quiver')
plot(pv(:,1),pv(:,2),'-*b')
h1 = streamline(xxx,yyy,u,v,pv(1,1)-1e4,pv(1,2));
h2 = streamline(xxx,yyy,u,v,pv(2,1),pv(2,2));
set(h1,'Color','none');
set(h2,'Color','none');
% plot(h1.XData(1:100:2700),h1.YData(1:100:2700),'*-','color',rgb('Forest Green'))
% plot(h2.XData(1:100:1300),h2.YData(1:100:1300),'*-','color',rgb('Forest Green'))
% pv_new = [h1.XData(1:100:2700)',h1.YData(1:100:2700)';flipud(h2.XData(1:100:1300)'),flipud(h2.YData(1:100:1300)');h1.XData(1), h1.YData(1)];
% plot(pv_new(:,1),pv_new(:,2))
% h = streamline(xi,yi,-u,-v,pv(:,1),pv(:,2));
% set(h,'Color',rgb('purple'));
% h = streamline(xi,yi,-u,-v,-2.95e5:.5e3:-2.9e5,-5.15e5:.5e3:-5.1e5);
% set(h,'Color',rgb('lavender'));
axis equal
save('pv_new','pv_new')

%%


thin_m = 0;
str = 'ISSM';
mapFile = 'gridFlowRiseC035.mat';
initializeInputs();
initializeModel();
tau_c = defineTau(str);
buildSystem();

[u,v] = measures_interp('velocity',xy(:,1),xy(:,2));
u = u/3.154e7;
v = v/3.154e7;

figure(2)
clf
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'tight');

nexttile(1)
contour(xxx,yyy,spd2, [10, 10] , 'k:','HandleVisibility','off')
hold on
contour(xxx,yyy,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xxx,yyy,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xxx,yyy,spd2, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)

outA = zeros(size(A*u));
outA(se_bound_c) = A(se_bound_c,:)*u;
trisurf(t_c,xy_c(:,1),xy_c(:,2),zeros(size(xy_c(:,1))),outA,...
               'edgecolor','none')
caxis([-1*max(abs(A*u)),max(abs(A*u))]);
colorbar
view(2)
colormap redblue

nexttile(2);
contour(xxx,yyy,spd2, [10, 10] , 'k:','HandleVisibility','off')
hold on
contour(xxx,yyy,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xxx,yyy,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xxx,yyy,spd2, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)

outB = zeros(size(A*u));
outB(se_bound_c) = B(se_bound_c,:)*v;
trisurf(t_c,xy_c(:,1),xy_c(:,2),zeros(size(xy_c(:,1))),outB,...
               'edgecolor','none')
caxis([-1*max(abs(B*v)),max(abs(B*v))]);
colorbar
view(2)
colormap redblue




