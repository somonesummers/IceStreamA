clear
if(~ismac)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
    run /home/groups/jsuckale/psummers/MATLAB/startup.m
end

load gridSiple5000;
figure(1)

[xy,t]=distmesh2d(@dpoly,@huniform,.01,[xmin,ymin;xmax,ymax]/1e5,pv/1e5,pv/1e5);
xy = xy*1e5;
save("gridSiple1000.mat");
