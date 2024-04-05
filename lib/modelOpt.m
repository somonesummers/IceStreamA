function [resid] = modelOpt(x0,str,fg1,fg2)
% Mapview code, bringing in Tau constraints from SedMap estimates.
% Requires download and install the following software: 
% Distmesh https://popersson.github.io/distmesh/index.html
% CVX http://cvxr.com/cvx/download/
% MEaSUREs matlab plug in https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
% BedMachine matlab plug in https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine

%% Initialization

mapFile = "gridFlowRiseA02.mat";
runType = 4;
% Load input files
thin_m = 0;
initializeInputs();

% Build model bed, geometery, functions
initializeModel();

%% Define \tau
% Define map of basal strength according to named scenario. 
tau_c = defineTau(str,x0);

%% Build System
buildSystem();

%% Thermomechanical coupling loop
for t_i = 1:100  
    % Thermocouple fields to update everyloop
    % Strain rate [s^-1]
        ep_dot = calcTrigridStrain(u,v,xy,dx); %returns intperolation object
        lambda  = calcAdvection(T,u,v,xy,dx/4,rho,C_p); %TODO better derivatives, analytic?

        % Brinkman number [ ]
        Br =@(x,y) 2*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y))).*((subplus(ep_dot(x,y)).^(nn+1))/A_m).^(1/nn);

        % Peclet number  [ ]
        Pe =@(x,y) rho*C_p.*Acc(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y))./K;

        % Vertical Peclet number  [ ]
        %La =@(x,y) lambda(x,y).*subplus(h_s_init(x,y)-h_b_init(x,y)).^2./(K*(T_m-T_s(x,y)));
        La =@(x,y) zeros(size(x));

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
    cap = 30^(-1/nn); %stability cap on enhancement
    e_new = (E_t(xy_c(:,1),xy_c(:,2))).^(-1/nn);
    e_new(e_new < cap) = cap;  % max enhancement is a min viscosity
    
    if(t_i == 1)
        enhance = e_new;
    end

    if(t_i ~= 1) %first step we don't relax, we use E = 1 everywhere (zero strain is also an options)
        enhance = (1-nu) * enhance + nu*e_new;
        res = norm(e_new - enhance) / norm(e_new);
%         disp("Residual: " + res); %Don't comment on this for opt search
        
        if (res < 1e-3) %check for thermal stabilization
            break; 
        end
    end

    speedUp = 1;
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
            u(ne_bound) == spd_BC_u_ne./3.154E7*speedUp;
            v(ne_bound) == spd_BC_v_ne./3.154E7*speedUp;
            u(sw_bound) == spd_BC_u_sw./3.154E7*speedUp;
            v(sw_bound) == spd_BC_v_sw./3.154E7*speedUp;
            u(nw_bound) == spd_BC_u_nw./3.154E7*speedUp;
            v(nw_bound) == spd_BC_v_nw./3.154E7*speedUp;
        minimize(obj)
    cvx_end
    if(~strcmp(cvx_status,"Solved"))
        disp("CVX has issues: " + cvx_status);
        if(~contains(cvx_status,"Solved"))
            break;
        end
    end
    % u and v are [m/s]    
    %% Visualization in loop 
    %(uncomment to see avg temp, enhancement, and Pe, Lambda, Br every loop
%     inLoopPlotting;
end

%% Vis out of loop
% spd2 = measures_interp('speed',xy(:,1),xy(:,2));
[u2,v2] = measures_interp('velocity',xy(:,1),xy(:,2)); %[m/yr]

if(sum(isnan(u2)) + sum(isnan(v2)) > 0)
    uFill = scatteredInterpolant(xy(~isnan(u2),1),xy(~isnan(u2),2),u2(~isnan(u2)));
    vFill = scatteredInterpolant(xy(~isnan(v2),1),xy(~isnan(v2),2),v2(~isnan(v2)));
    u2(isnan(u2)) = uFill(xy(isnan(u2),1),xy(isnan(u2),2));
    v2(isnan(v2)) = vFill(xy(isnan(v2),1),xy(isnan(v2),2));    
    clear uFill vFill
end

spd2 = sqrt(u2.^2 + v2.^2);
spd_star = sqrt(v.^2+u.^2)*3.154E7;
resid1 = sum((u*3.154e7-u2).^2);
resid2 = sum((v*3.154e7-v2).^2);
resid3 = sum(log((spd_star+1)./(spd2+1)).^2);
% disp((resid1+resid2)./resid3);
gamma = 5e3; % 5e3 is even weight roughly
resid_xy = ((u*3.154e7-u2).^2 + (v*3.154e7-v2).^2 + gamma*log((spd_star+1)./(spd2+1)).^2).*F'/1e7; %% weight each by F matrix
resid_xy(xy(:,1)>-3e5) = 0; %don't care about far right region of tongue.
resid_xy(xy(:,2)>-4.5e5) = 0; %don't care about upper region.
resid = sum(resid_xy);

if exist('fg1','var') == 1
    set(0,'CurrentFigure',fg1);
else 
    fg1 = figure(1);
end
clf
trisurf(t,xy(:,1),xy(:,2),resid_xy,...
       'edgecolor','none')
title('Residual Contribution')
xlabel('X')
ylabel('Y')
colorbar
view(2)
% axis equal
drawnow
% figure
% clf
% trisurf(t_c,xy_c(:,1),xy_c(:,2),enhance.^(-nn),...
%        'edgecolor','none')
% title('Enhancement')
% view(2)
% axis equal
% xlabel('X')
% ylabel('Y')
% colorbar

% figure
% clf
% trisurf(t,xy(:,1),xy(:,2),zeros(size(xy(:,1))),phi_init(xy(:,1),xy(:,2)),...
%        'edgecolor','none','facecolor','interp')
% title('\phi')
% view(2)
% axis equal
% hold on
% tricontour(t,xy(:,1),xy(:,2),spd2,[10,30,100,300,1000,3000]);
% xlabel('X')
% ylabel('Y')
% colorbar

% figure('Position', [0 0 1200 600]);
% clf
% 
% subplot(121)
% ep = calcTrigridStrain(u2/3.154E7,v2/3.154E7,xy,dx);
% trisurf(t_c,xy_c(:,1),xy_c(:,2),log10(subplus(ep(xy_c(:,1),xy_c(:,2)))),...
%        'edgecolor','none')
% colorbar
% view(2)
% title('Observed Strain')
% 
% subplot(122)
% ep = calcTrigridStrain(u,v,xy,dx);
% trisurf(t_c,xy_c(:,1),xy_c(:,2),log10(subplus(ep(xy_c(:,1),xy_c(:,2)))),...
%        'edgecolor','none')
% view(2)
% colorbar
% title('Model Strain')



%figure('Position', [0 0 1200 600]);
if exist('fg2','var') == 1
    set(0,'CurrentFigure',fg2);
else 
    fg2 = figure(3);
end
clf
set(gcf, 'defaultAxesFontName', 'National Park Regular',...
'defaultTextFontName', 'National Park Regular');
sgtitle(str + " X_0 = " + x0 +"; Time " + paulTime());
subplot(141)

trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),log10(spd2),...
       'edgecolor','none')
hold on
trisurf(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),...
       'edgecolor','black','facecolor','none')

title('Speed of Measures')
xlabel('X')
ylabel('Y')
caxis([1 2.6]);
colorbar
view(2)
% axis equal


subplot(142)
trisurf(t,xy(:,1),xy(:,2),h_s_init(xy(:,1),xy(:,2)),log10(sqrt(u.^2 + v.^2)*3.154E7),...
       'edgecolor','none')
caxis([1 2.6]);   
title('Speed')
xlabel('X')
ylabel('Y')
colorbar
view(2)
% axis equal

subplot(143)
trisurf(t,xy(:,1),xy(:,2),tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2),...
       'edgecolor','none')
hold on
trisurf(t,xy(:,1),xy(:,2),h_b_init(xy(:,1),xy(:,2)),...
       'edgecolor','black','facecolor','none')
colorbar
caxis([00e3 150e3]);
colormap(gca, Cmap/255.0)
title('Basal \tau')
xlabel('X')
ylabel('Y')
view(2)
% axis equal

subplot(144)
trisurf(t_c,xy_c(:,1),xy_c(:,2),df,...
    'edgecolor','none');
    %'facecolor','interp');
title('Driving force')
xlabel('X')
ylabel('Y')
colorbar
caxis([0e3 150e3]);
colormap(gca, Cmap/255.0)
view(2)
% axis equal

end
