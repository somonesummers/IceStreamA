clear; clc; close all

if(~ismac)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
    run /home/groups/jsuckale/psummers/MATLAB/startup.m
end
fID = fopen('log.txt','w');

x0 = [1.4,.5];
str = "ISSM Shift";
fg1 = figure(1);
fg2 = figure(2);
fun = @(x)modelOpt(x,str,fg1,fg2);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e-2    ,'TolFun',1e4);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

save("optOutput_"+ str); 
