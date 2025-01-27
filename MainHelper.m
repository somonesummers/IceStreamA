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
addpath('grids','lib')
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
thinToRun = [50,40,30]; 
speedUpToRun = [1];
NToRun = [1]; % 0 = constant Eff, 1 = adjust Eff with HAF 
runType = 5; %1 = dH explicit, 2 = DhDt based, 3 =  Golledge based, 4 = No Lakes, 5 mass bal
%N_adjust = 1; %1 = change Tau with N, 0 Tau fixed
for j = 1:length(NToRun)
    for i = 1:length(thinToRun)
        clearvars -except nameToRun mapsToRun thinToRun speedUpToRun i j fID runType N_adjust NToRun
        thin_m = thinToRun(i);
        N_adjust = NToRun(j);
        speedUp = speedUpToRun(1);
        str = nameToRun(1);
        mapFile = mapsToRun(1);
        ModelRunner;
    end
end
finalString = 'exit from MainHelper successfully ' + paulTime();
disp(finalString);
fprintf(fID,finalString);
