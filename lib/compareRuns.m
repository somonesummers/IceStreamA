clear; close all

%% Load data

for i = 4
    if(i == 1)
    data1 = load('data/data_gridInstitute24000MixedElmernoAdvect.mat');
    data2 = load('data/data_gridInstitute24000MixedElmer_streamnoAdvect.mat');
    elseif(i ==2)
    data1 = load('data/data_gridInstitute24000ISSM_centernoAdvect.mat');
    data2 = load('data/data_gridInstitute24000ISSM_center_streamnoAdvect.mat');
    elseif(i ==3)
    data1 = load('data/data_gridInstitute24000ELMER_centernoAdvect.mat');
    data2 = load('data/data_gridInstitute24000ELMER_center_streamnoAdvect.mat');
    else
    data1 = load('data/data_gridSipleSmall7500ISSMbedmap50.mat');
    data2 = load('data/data_gridSipleSmall7500ISSMbedmap0.mat');
    end
    runPlotting(data1,data2);
end
    


function [] = runPlotting(data1,data2)    


ftsize = 18;

%% Plotting speeds
figure
clf
subplot(131)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        data1.h_s,(sqrt(data1.u.^2 + data1.v.^2)*3.154E7...
        - sqrt(data2.u.^2 + data2.v.^2)*3.154E7)...
           ,'edgecolor','none');
    title('Speed diff (1-2)')
    xlabel('X')
    ylabel('Y')
    colorbar
    view(2)
    axis equal
    

subplot(132)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        data1.h_s,log10(sqrt(data1.u.^2 + data1.v.^2)*3.154E7)...
           ,'edgecolor','none');
    title('Speed 1')
    xlabel('X')
    ylabel('Y')
    colorbar

    caxis([1 2.8]);
    view(2)
    axis equal

subplot(133)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        data1.h_s, log10(sqrt(data2.u.^2 + data2.v.^2)*3.154E7)...
           ,'edgecolor','none');
    title('Speed 2')
    xlabel('X')
    ylabel('Y')
    colorbar
    caxis([1 2.8]);
    view(2)
    axis equal

%% plotting Tau
figure
clf
subplot(141)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        data1.tau_c(data1.xy(:,1),data1.xy(:,2),data1.u,data1.v)./norms([data1.u,data1.v],2,2)/1e3...
        - data2.tau_c(data2.xy(:,1),data2.xy(:,2),data2.u,data2.v)./norms([data2.u,data2.v],2,2)/1e3...
           ,'edgecolor','none');
    title('\tau diff (1-2)')
    xlabel('X')
    ylabel('Y')
    colorbar
    view(2)
    axis equal

subplot(142)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        data1.h_s,data1.tau_c(data1.xy(:,1),data1.xy(:,2),data1.u,data1.v)./norms([data1.u,data1.v],2,2)/1e3...
           ,'edgecolor','none');
    title('\tau 1')
    xlabel('X')
    ylabel('Y')
    colorbar
    caxis([0 100]);
    view(2)
    axis equal

subplot(143)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        data1.h_s, data2.tau_c(data2.xy(:,1),data2.xy(:,2),data2.u,data2.v)./norms([data2.u,data2.v],2,2)/1e3...
           ,'edgecolor','none');
    title('\tau 2')
    xlabel('X')
    ylabel('Y')
    colorbar
    caxis([0 100]);
    view(2)
    axis equal
subplot(144)
    trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        (data1.tau_c(data1.xy(:,1),data1.xy(:,2),data1.u,data1.v)./norms([data1.u,data1.v],2,2)...
        - data2.tau_c(data2.xy(:,1),data2.xy(:,2),data2.u,data2.v)./norms([data2.u,data2.v],2,2))...
        ./(data1.tau_c(data1.xy(:,1),data1.xy(:,2),data1.u,data1.v)./norms([data1.u,data1.v],2,2)),'edgecolor','none');
    title('\tau diff (1-2) percent')
    xlabel('X')
    ylabel('Y')
    colorbar
%     caxis([-.25 .25])
    view(2)
    axis equal
%% Just speed diff
figure
clf
trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),...
        0*data1.h_s,sqrt(data1.u.^2 + data1.v.^2)*3.154E7...
        - sqrt(data2.u.^2 + data2.v.^2)*3.154E7...
           ,'edgecolor','none','facecolor','interp');
hold on
% Data 1 Contour
% [C, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
%     sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[100 300 3000]);
% for j = 1:numel(H)
% H(j).EdgeColor = rgb('black');
% end
[C, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
    sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[30 30]);
for j = 1:numel(H)
H(j).EdgeColor = rgb('black');
H(j).LineStyle = '--';
H(j).LineWidth = 3;
end
% [C, H] = tricontour(data1.t,data1.xy(:,1),data1.xy(:,2),...
%     sqrt(data1.u.^2 + data1.v.^2)*3.154E7,[1000 1000]);
% for j = 1:numel(H)
% H(j).EdgeColor = rgb('black');
% H(j).LineWidth = 2;
% end

% Data 2 Contour
% [C, H] = tricontour(data2.t,data2.xy(:,1),data2.xy(:,2),...
%     sqrt(data2.u.^2 + data2.v.^2)*3.154E7,[100 300 3000]);
% for j = 1:numel(H)
% H(j).EdgeColor = rgb('grey');
% end
[C, H] = tricontour(data2.t,data2.xy(:,1),data2.xy(:,2),...
    sqrt(data2.u.^2 + data2.v.^2)*3.154E7,[30 30]);
for j = 1:numel(H)
H(j).EdgeColor = rgb('grey');
H(j).LineStyle = '--';
H(j).LineWidth = 3;
end
% [C, H] = tricontour(data2.t,data2.xy(:,1),data2.xy(:,2),...
%     sqrt(data2.u.^2 + data2.v.^2)*3.154E7,[1000 1000]);
% for j = 1:numel(H)
% H(j).EdgeColor = rgb('grey');
% H(j).LineWidth = 2;
% end


CT2 = flipud(cbrewer('seq','Reds',256));
load fire.mat
CT = 1-fire/max(max(fire));
title('Speed diff (1-2) ' + data1.str)
% colormap(CT)
colormap redblue
xlabel('Easting')
ylabel('Northing')
c = colorbar;
% caxis([-225 0])
caxis([-30 30])
c.Label.String = '\Delta Ice Speed [m/yr]';
view(2)
f = gca;
f.XAxis.FontSize = ftsize;
f.YAxis.FontSize = ftsize;
c.FontSize = ftsize;
view(2)
axis equal

% saveVect("figs/Diff" + data1.str + "redblue")

%% Compare to anti-flow model
figure
subplot(121)
trisurf(data1.t,data1.xy(:,1),data1.xy(:,2),zeros(size(data1.xy(:,1))),(sqrt(data1.u.^2 + data1.v.^2)*3.154E7),...
       'edgecolor','none')   
caxis([10  675])
hold on
% quiver(data1.xy(:,1),data1.xy(:,2),data1.u,data1.v)
title('Speed')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

subplot(122)
trisurf(data2.t,data2.xy(:,1),data2.xy(:,2),zeros(size(data2.xy(:,1))),(sqrt(data2.u.^2 + data2.v.^2)*3.154E7),...
       'edgecolor','none')   
caxis([10 675])
hold on
% quiver(data2.xy(:,1),data2.xy(:,2),data2.u,data2.v)
title('Speed')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

end
