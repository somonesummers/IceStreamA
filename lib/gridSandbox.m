% File used to test/sandbox out functions on the grid. Mostly for work importing from external datasets and make plots
%% Make Grid
clear
% close all
addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
run /home/groups/jsuckale/psummers/MATLAB/startup.m

dx = 1e3;    %dx: nominal grid spacing [m]

fname = ('ice-stream-a/ice-stream-a-domain.geojson');
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
clear str
xbox = val.features.geometry.coordinates(:,:,1);
ybox = val.features.geometry.coordinates(:,:,2);
%     small
% xbox = [-6.2457 -4.0000   -2.5000   -4.0000  -6.2457]*1e5;
% ybox = [-4.600 -2.7000   -4.9350   -6.3000   -4.600]*1e5;
%    xsm
%    xbox = [-4.517 -3.260   -2.4500   -3.742   -4.517]*1e5;
%    ybox = [-4.000 -3.658   -5.0350   -5.600   -4.000]*1e5;
xmax =  max(xbox);
xmin = min(xbox);
ymax =  max(ybox);
ymin =  min(ybox);


%     theta = pi*.20; %[rad] angle of rotation
% Mesh Generation
clf
pv = [xbox; ybox]';
%     pv = makePerimeter(pv,dx); %Places points on the perimeter for distMesh to use
figure(2)
[xy,t] = distmesh2d(@dpoly,@huniform,dx/1e5,[xmin,ymin;xmax,ymax]/1e5,pv/1e5,pv/1e5);
% scale to real grid size
xy = xy*1e5;
%grab boundry points
se_bound = (ybox(4)-ybox(3))/(xbox(4)-xbox(3))*xy(:,1) - xy(:,2)  > (ybox(4)-ybox(3))/(xbox(4)-xbox(3))*xbox(3) - ybox(3)-dx/3;
ne_bound = (ybox(3)-ybox(2))/(xbox(3)-xbox(2))*xy(:,1) - xy(:,2)  < (ybox(3)-ybox(2))/(xbox(3)-xbox(2))*xbox(2) - ybox(2)+dx/3;
nw_bound = (ybox(2)-ybox(1))/(xbox(2)-xbox(1))*xy(:,1) - xy(:,2)  < (ybox(2)-ybox(1))/(xbox(2)-xbox(1))*xbox(1) - ybox(1)+dx/3;
sw_bound = (ybox(1)-ybox(4))/(xbox(1)-xbox(4))*xy(:,1) - xy(:,2)  > (ybox(1)-ybox(4))/(xbox(1)-xbox(4))*xbox(4) - ybox(4)-dx/3;


save("gridSiple" + dx + ".mat");

%%
if(ismac)
    figure
    spd = measures_interp('speed',xy(:,1),xy(:,2));
    trisurf(t,xy(:,1),xy(:,2),zeros(size(spd)),(spd),...
           'edgecolor','none')
    view(2)
    hold on
    plot(pv(:,1),pv(:,2),'r-','linewidth',2)
end



