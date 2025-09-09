
clear 
% baseFile = "data/datagridFlowRiseA02UpBCISSMMb_DhDt0SpeedUp0.mat";
baseFile = "data/spdChange/data_NgridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
% groupName = 'fluxBySpeed';
% cases = [0,1,2,3,4,5];
groupName = 'fluxByThinning';
cases = [0,10,20,30,40,50];

flux = zeros(numel(cases),1);

for j = 1:numel(cases)
%     baseFile = "data/datagridFlowRiseA02ISSMNoLakes_DhDt0SpeedUp0.mat";
    newFile = strrep(baseFile,"Dt0","Dt" + cases(j));
    % newFile = strrep(baseFile,"Up0","Up" + cases(j));

    data = load(newFile);
    flux(j) = mean(data.se_bound.*(data.u.^2 + data.v.^2).^(.5).*data.h)*3.154e7;
    clear data
end

figure(1)
plot(cases*.3,flux,'-o')
xlabel('Thickness change [m]')
% plot(cases*10,flux,'-o')
% xlabel('Speed Up [%]')
ylabel('unit Flux out bottom [m^2/yr]')

setFontSize(18)

fig = gcf;
% labelTiledLayout(fig, 18)
savePng("figs/paper/" + groupName);