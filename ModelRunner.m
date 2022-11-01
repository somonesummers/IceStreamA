% Requires download and install the following software: 
% Distmesh https://popersson.github.io/distmesh/index.html
% CVX http://cvxr.com/cvx/download/
% MEaSUREs matlab plug in + data files
% https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
% BedMachine matlab plug in + data files
% https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine

% clc
% clear
% close all

% fID = fopen('log.txt','w');

%% Initialization
% Scenario to run if running one at a time comment out below, run this file
% directly
saveData = true;
% str = 'Mixed';
% mapFile = 'gridInstitute5000.mat';

% Comment so we know what's happening, thats always nice.
disp("Running " + str + " " + thin_m + " now...");
fprintf(fID,'\tRunning %s %d now...\n',str,thin_m);
if(~saveData)
    disp("Data will not be saved");
    fprintf(fID,'Data will not be saved');
end

% Load input files
initializeInputs();

% Build model bed, geometery, functions
initializeModel();

%% Define \tau
% Define map of basal strength according to named scenario. 
tau_c = defineTau(str);
% old_tau = load('data/data_gridSipleXXSmall5000ISSM Shift2bedmap0.mat');
% tau_c = old_tau.tau_c;
%% Build System
buildSystem();


%% Thermomechanical coupling loop
for t_i = 1:100  
    % Thermocouple fields to update everyloop
    % Strain rate [s^-1]
        if t_i == 1
            [u,v] = measures_interp('velocity',xy(:,1),xy(:,2));
            u = u/3.514e7;
            v = v/3.514e7;
            if(sum(isnan(u)) + sum(isnan(v)) > 0)
                uFill = scatteredInterpolant(xy(~isnan(u),1),xy(~isnan(u),2),u(~isnan(u)));
                vFill = scatteredInterpolant(xy(~isnan(v),1),xy(~isnan(v),2),v(~isnan(v)));
                u(isnan(u)) = uFill(xy(isnan(u),1),xy(isnan(u),2));
                v(isnan(v)) = vFill(xy(isnan(v),1),xy(isnan(v),2));
                clear uFill vFill
            end
        end
        ep_dot = calcTrigridStrain(u,v,xy,dx); %returns intperolation object
        
        if(true)
            T_calc = T;
        else
            T_calc = T_bar(xy(:,1),xy(:,2));
        end
        lambda  = calcAdvection(T_calc,u,v,xy,dx,rho,C_p); 

        % Brinkman number [ ]
        Br =@(x,y) 2*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y))).*((subplus(ep_dot(x,y)).^(nn+1))/A_m).^(1/nn);

        % Peclet number  [ ]
        Pe =@(x,y) rho*C_p.*Acc(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y))./K;

        % Horizontal Peclet number  [ ]
%         La =@(x,y) lambda(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y)));
        La =@(x,y) zeros(size(x)); %exclude advection
        
        % Critical Strain [s^-1]
        ep_star =@(x,y) ((La(x,y)/2 + ((Pe(x,y).^2)/2)./(Pe(x,y)-1+exp(-Pe(x,y))))).^(nn/(nn+1))...
        .*(K*(T_m-T_s(x,y))./(A_m.^(-1/nn).*(subplus(h_s_init(x,y)-h_b_init(x,y))).^2)).^(nn/(nn+1));

        % Temp profile at xy [K]
        t_z =@(x,y) tempProfile(ep_dot(x,y),ep_star(x,y),Pe(x,y),Br(x,y),La(x,y),T_s(x,y),T_m,dz); 

        % Enhancement Factor []
        E_t =@(x,y) depthIntEnhancement(t_z(x,y),a.^(-3),dz);

        % Mean Temp [K]
        T_bar = @(x,y) trapz(t_z(x,y),2)*dz;

    % Calc Enhancement Factors, relax into solution. Have max value for
    % stabilization
    cap = 20^(-1/nn); %stability cap on enhancement
    e_new = (E_t(xy_c(:,1),xy_c(:,2))).^(-1/nn);
%     e_new(e_new < cap) = cap;  % max enhancement is a min viscosity
    
    if(t_i == 1)
        enhance = e_new;
    end

    if(t_i ~= 1) %first step we don't relax, we use E = 1 everywhere (zero strain is also an options)
        enhance = (1-nu) * enhance + nu*e_new;
        res = norm(e_new - enhance) / norm(e_new);
        disp("Residual: " + res);
        fprintf(fID,'\t\t[%d]Residual: %f \n',t_i,res);
        if (res < 1e-3) %check for thermal stabilization
            break; 
        end
    end
    %% Visualization in loop 
    %(uncomment to see avg temp, enhancement, and Pe, Lambda, Br every loop
    if(ismac)
        inLoopPlotting;
    end
    %% Solve
    % Unused BCs
    %       v(xy(:,2) > ymax - dx/2) == 0;
    %       v(xy(:,2) < ymin - dx/2) == 0;
    
    % Options include     cvx_precision low, cvx_begin quiet
    % CVX may throw a warning about non-empty problems here, that is OK.
    cvx_begin quiet
        variables u(nN) v(nN)
        obj = 2.*a./p.*sum(enhance.*h_av.*tau_area.*pow_pos(norms([A*u,B*v,1/2*(B*u+A*v)],2,2),p)) + ...
              F*tau_c(xy(:,1),xy(:,2),u,v) + ...
              rho*g*sum(h_av.*((A*h_s).*(D*u) + (B*h_s).*(D*v)));
        subject to
            u(se_bound) == spd_BC_u_se./3.154E7;
            v(se_bound) == spd_BC_v_se./3.154E7;
            u(nw_bound) == spd_BC_u_nw./3.154E7;
            v(nw_bound) == spd_BC_v_nw./3.154E7;
            u(ne_bound) == spd_BC_u_ne./3.154E7;
            v(ne_bound) == spd_BC_v_ne./3.154E7;
            u(sw_bound) == spd_BC_u_sw./3.154E7;
            v(sw_bound) == spd_BC_v_sw./3.154E7;
            
        minimize(obj)
    cvx_end
    if(~strcmp(cvx_status,"Solved"))
        disp("CVX has issues: " + cvx_status);
        fprintf(fID,'\t\tCVX has issues: %s\n',cvx_status);
        if(~contains(cvx_status,"Solved"))
            break;
        end
    end
    % u and v are [m/s]    
    if(t_i == 1)
        u_init = u;
        v_init = v;
    end
    
end
if(ismac)
    clear fg1 fg2
end
%% Save data to data file
mpClean = erase(mapFile, [".mat","workingGrid_"]);
if(saveData)
    save("data/data_" + mpClean + str + "bedmap_2_" + thin_m + ".mat");
else
    warning('Data not being saved');
    Disp('Data not being saved');
end

%% Vis out of loop
spd2 = measures_interp('speed',xy(:,1),xy(:,2)); %[m/yr]

figure('Position', [0 0 1200 600]);
clf
sgtitle(str);
subplot(141)

trisurf(t,xy(:,1),xy(:,2),zeros(size(spd2)),(spd2),...
       'edgecolor','none')
hold on
title('Speed of Measures')
xlabel('X')
ylabel('Y')
caxis([10  675])
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
view(2)
axis equal

subplot(142)
trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),(sqrt(u.^2 + v.^2)*3.154E7),...
       'edgecolor','none')   
caxis([10  675])
title('Speed')
xlabel('X')
ylabel('Y')
colorbar
f = gca;
f.ColorScale = 'log';
view(2)
axis equal

subplot(143)
trisurf(t,xy(:,1),xy(:,2),tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2),...
       'edgecolor','none')
% hold on
% trisurf(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),...
%        'edgecolor','black','facecolor','none')
colorbar
caxis([0e3 150e3]);
colormap(gca, Cmap/255.0)
title('Basal \tau')
xlabel('X')
ylabel('Y')
view(2)
axis equal

subplot(144)
trisurf(t_c,xy_c(:,1),xy_c(:,2),df,...
    'edgecolor','none');
title('Driving force')
xlabel('X')
ylabel('Y')
colorbar
caxis([0e3 150e3]);
colormap(gca, Cmap/255.0)
view(2)
axis equal

int_x = scatteredInterpolant(xy(:,1),xy(:,2),u,'natural','none');
int_y = scatteredInterpolant(xy(:,1),xy(:,2),v,'natural','none');
[msr_x, msr_y] = measures_interp('velocity',Xi,Yi);

spdModel = (sqrt(u.^2 + v.^2)*3.154E7);

figure()
quiver(Xi,Yi,msr_x,msr_y)
hold on 
quiver(Xi,Yi,int_x(Xi,Yi),int_y(Xi,Yi))
scatter(xy(spdModel > spd2*2,1),xy(spdModel > spd2*2,2),'r','filled')

if(exist('u_init','var'))
    figure
    trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),((sqrt(u.^2 + v.^2))-(sqrt(u_init.^2 + v_init.^2)))*3.154E7,...
           'edgecolor','none')
    title('Coupling Speedup')
    xlabel('X')
    ylabel('Y')
    colorbar
    colormap redblue
    caxis([-300 300])
    view(2)
    axis equal
end
