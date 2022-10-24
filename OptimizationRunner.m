clear; clc; close all

if(~ismac)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
    run /home/groups/jsuckale/psummers/MATLAB/startup.m
end
fID = fopen('log.txt','w');

x0 = [1.4,.7];
str = "ISSM Shift";

fun = @(x)modelOpt(x,str);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e3    ,'TolFun',1e6);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

save("optOutput_"+ str); 
