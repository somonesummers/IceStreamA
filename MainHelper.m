%% Main file to run, allows for multiple scenarios to be run in series
% Requires download and install the following: 
% 
% Distmesh 
%    https://popersson.github.io/distmesh/index.html
% CVX 
%    http://cvxr.com/cvx/download/
% MEaSUREs matlab plug in + data files
%   https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
% BedMachine matlab plug in + data files
%   https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine


clear
close all
clc
addpath('grids')
if(~ismac)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
    run /home/groups/jsuckale/psummers/MATLAB/startup.m
end
fID = fopen('log.txt','w');
% Name of scenarios to run, only 1 map file used here.
nameToRun = ["ISSM"];
mapsToRun = ["gridFlowRiseA02.mat"];
thinToRun = [0,10,20,30,40,50,60,70,80,100,-10,-20,-30];
speedUpToRun = [1];
for j = 1:length(speedUpToRun)
    for i = 1:length(thinToRun)
%         clearvars -except nameToRun mapsToRun thinToRun speedUpToRun i j fID
        thin_m = thinToRun(i);
        speedUp = speedUpToRun(j);
        str = nameToRun(1);
        mapFile = mapsToRun(1);
        ModelRunner;
    end
end
disp('exit from MainHelper successfully');
fprintf(fID,'exit from MainHelper successfully');