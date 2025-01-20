clear
if(~ismac)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
    run /home/groups/jsuckale/psummers/MATLAB/startup.m
end
% siple
% xbox = [-6.2457 -4.0000   -2.5000   -4.0000  -6.2457]*1e5;
% ybox = [-4.600 -2.7000   -4.9350   -6.3000   -4.600]*1e5;
% xsmall
% xbox = [-4.517 -3.260   -2.4500   -3.742   -4.517]*1e5;
% ybox = [-4.000 -3.658   -5.0350   -5.600   -4.000]*1e5;
% All of Ridge
xbox = [-4.9500  -3.6000   -2.7000   -4.1000   -4.9500]*1e5;
ybox = [-3.8000  -3.5300   -5.2250   -5.6700   -3.8000]*1e5;

xmax =  max(xbox);
xmin = min(xbox);
ymax =  max(ybox);
ymin =  min(ybox);

pv = [xbox; ybox]';

figure(1)

[xy,t]=distmesh2d(@dpoly,@huniform,.05,[xmin,ymin;xmax,ymax]/1e5,pv/1e5,pv/1e5);
xy = xy*1e5;
[u,v] = measures_interp('velocity',xy(:,1),xy(:,2));

dx = 1e3;
x = xmin:dx:xmax;
y = ymin:dx:ymax;
[xx,yy] = meshgrid(x,y);
spd = measures_interp('speed',xx,yy);
figure
surf(xx,yy,zeros(size(spd)),log10(spd),'edgecolor','none')
hold on 
view(2)
plot(pv(:,1),pv(:,2),'r-')
axis equal

[uu,vv] = measures_interp('velocity',xx,yy);


uu  = cleanNan(xx,yy, uu);
vv  = cleanNan(xx,yy, vv);

figure
quiver(xx,yy,uu./spd.*log10(spd),vv./spd.*log10(spd))

figure
quiver(xx,yy,uu./spd,vv./spd)
hold on
contour(xx,yy,spd,[10,30,100,300,1000],'k')

u   = cleanNan(xy(:,1),xy(:,2),  u);
v   = cleanNan(xy(:,1),xy(:,2),  v);

figure
trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),sqrt(u.^2+v.^2),...
       'edgecolor','none')  
colorbar
view(2)
%%
ftsize = 16;
figure('Position',[300 300 700 850])

uu = imgaussfilt(uu,2);
vv = imgaussfilt(vv,2);
alpha = atan(vv./uu);
[ux ,  uy] = gradient(uu,dx,dx);
[vx ,  vy] = gradient(vv,dx,dx);
e_xy = .5 * (uy + vx);

e_shr = (vy-ux).*cos(alpha).*sin(alpha) + e_xy.*(cos(alpha).^2 - (sin(alpha).^2));

e_shr(isnan(e_shr)) = 0;
surf(x,y,zeros(size(xx)),e_shr,...
     'edgecolor','none','facecolor','interp')
hold on
plot(pv(:,1),pv(:,2),'k-','LineWidth', 4)
contour(x,y,spd, [10, 10] , 'k:','HandleVisibility','off')
contour(x,y,spd, [30, 30] , 'k--','HandleVisibility','off')
contour(x,y,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(x,y,spd, [1000, 1000] , 'k-','HandleVisibility','off','LineWidth',2)
view(2)
CT = cmocean('curl');
colormap(CT)
axis equal
ylabel('Northing [m]')
xlabel('Easting [m]')
c = colorbar;
caxis([-0.005 0.005])
c.Label.String = 'Shear Strain Rate [yr^{-1}]';
c.FontSize = ftsize;   
view(2)
f = gca;
f.XAxis.FontSize = ftsize-2;
f.YAxis.FontSize = ftsize-2;
view(2)
axis equal
title('Observed Strain Rate','fontsize',ftsize)

%%
% ep  = calcTrigridStrain(u,v,xy,3e3);
    dx =2e3;
    xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
    yi = (min(xy(:,2))-dx:dx:max(xy(:,2))+dx)';
    [Xi,Yi] = ndgrid(xi,yi);
    us = scatteredInterpolant(xy(:,1),xy(:,2),u);
    vs = scatteredInterpolant(xy(:,1),xy(:,2),v);

    % grid to find strain rates
    ug = us(Xi,Yi); % [m/s]
    vg = vs(Xi,Yi); % [m/s]

    % Smooth
%     ug = imgaussfilt(ug,1);
%     vg = imgaussfilt(vg,1);
    
    % Calc Strain on grid
    [n,m] = size(ug);
    strain = zeros(n-1,m-1,2,2);
    strain(:,:,1,1) = diff(vg(1:end-1,:),1,2)/dx;
    strain(:,:,2,2) = diff(ug(:,1:end-1))/dx;
    strain(:,:,1,2) = .5 .* (diff(ug(1:end-1,:),1,2) + diff(vg(:,1:end-1)))/dx;
    strain(:,:,2,1) =  strain(:,:,1,2);
%   Total  Effective Strain
    EffStrain = sqrt(.5*(strain(:,:,1,1).^2 + strain(:,:,2,2).^2 +...
            (strain(:,:,1,1) + strain(:,:,2)).^2) + strain(:,:,1,2).^2);
    EffStrain(isnan(EffStrain)) = 1e-14; % this doesn't happen anymore Yay!
    EffStrain = imgaussfilt(EffStrain,5);
    ep = griddedInterpolant(Xi(2:end,2:end)-dx/2, Yi(2:end,2:end)-dx/2, EffStrain,'linear','nearest');




fun =@(x,y) min(subplus(log(ep(x*1e5,y*1e5))+9),5);
figure
scatter(xy(:,1),xy(:,2),[],fun(xy(:,1)/1e5,xy(:,2)/1e5),'filled')
colorbar


% figure(2)
% scatter(xx(:),yy(:),[],lgSpd(:),'filled')
% colorbar

% %%
% fd=@(p) dpoly(p,pv/1e5);
% % fh=@(p) 10-log(measures_interp('speed',p(:,1)*1e5,p(:,2)*1e5)+1);
% fh=@(p) ones(size(p,1),1) - fun(p(:,1),p(:,2))/6;
% 
% figure
% edgeLength = .03;
% [xy,t]=distmesh2d(fd,fh,edgeLength,[xmin,ymin;xmax,ymax]/1e5,pv/1e5);
% 
% xy = xy*1e5;
% se_bound = (ybox(4)-ybox(3))/(xbox(4)-xbox(3))*xy(:,1) - xy(:,2)  > (ybox(4)-ybox(3))/(xbox(4)-xbox(3))*xbox(3) - ybox(3)-dx/3;
% ne_bound = (ybox(3)-ybox(2))/(xbox(3)-xbox(2))*xy(:,1) - xy(:,2)  < (ybox(3)-ybox(2))/(xbox(3)-xbox(2))*xbox(2) - ybox(2)+dx/3;
% nw_bound = (ybox(2)-ybox(1))/(xbox(2)-xbox(1))*xy(:,1) - xy(:,2)  < (ybox(2)-ybox(1))/(xbox(2)-xbox(1))*xbox(1) - ybox(1)+dx/3;
% sw_bound = (ybox(1)-ybox(4))/(xbox(1)-xbox(4))*xy(:,1) - xy(:,2)  > (ybox(1)-ybox(4))/(xbox(1)-xbox(4))*xbox(4) - ybox(4)-dx/3;
% 
% figure
% clear u
% u = measures_interp('speed',xy(:,1),xy(:,2));
% uFill = scatteredInterpolant(xy(~isnan(u),1),xy(~isnan(u),2),u(~isnan(u)));
% u(isnan(u)) = uFill(xy(isnan(u),1),xy(isnan(u),1));
%    
% trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),u,...
%        'edgecolor','none')   
% view(2)
% colorbar   

% save("gridRefinedXSM" + strrep(string(edgeLength),"0.","") + ".mat");
