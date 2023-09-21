clear; clc; close all

if(~ismac)
    %Start-up business on sherlock is hard 
    addpath('lib') 
    addpath(genpath('/home/groups/jsuckale/psummers/MATLAB'))
    run /home/groups/jsuckale/psummers/MATLAB/startup.m
end
fID = fopen('log.txt','w');

x0 = [58e3];
str = "Uniform";
fg1 = figure(1);
fg2 = figure(2);
fun = @(x)modelOpt(x,str,fg1,fg2);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e-1    ,'TolFun',1e4);	
[x,fval,exitflag,output] = fminsearch(fun,x0,options) 

clear fg1 fg2
save("optOutput_"+ str); 
