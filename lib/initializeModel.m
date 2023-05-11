% Initialize Model parameters, must be done after loading inputs

%% Model Vars
nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes
e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges
res = 1; %start residual arb high number

%% Physical parameters
a     = 3.5e-25^(-1/3);     % a:     flow parameter pre-factor [Pa s^1/3] @-10C from cuffey 
nn    = 3;                  % Glens law power
p     = 4/3;                % p:     flow parameter power [ ]
g     = 9.81;               % g:     acceleration due to gravity [m/s^2]
nu    = .5;                 % Thermal relaxation parameter [ ]
rho   = 917;                % rho:   density of ice [kg/m^3]
rho_w = 1000;               % rho_w: density of water [kg/m^3]
C_p   = 2050;               % specific heat of ice [J/Kg/K]
K     = 2.1;                % thermal conductivity of ice [W/m/K]
A_m   = 2.4e-24;            % Meyer's prefactor [Pa^-3 s^-1]
T_m   = 273;                % Ice melting point [k] 
enhance = ones(nE,1);       % Thermal enhancement factor [ ]
E_man = 1; %(1/6)^(-1/3);       % Manual Adjustment to thermocoupling enhancement

% initialize to zero velocity case [m/s] 
u = zeros(size(xy(:,1))); 
v = u;

dz = .1;  %vertical resolution of thermal depth profiles (frac of H) [ ]

%% Import bedMachine data and smooth 
% import data on higher resolution square grid that is larger than the model grid by
% 'overgrab'
xmax =  max(pv(:,1));
xmin = min(pv(:,1));
ymax =  max(pv(:,2));
ymin =  min(pv(:,2));

overgrab = 20;
xi = xmin-dx*overgrab:dx/2:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx/2:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% % Raw fields
% % bm_b =  bedmap2_interp(Xi,Yi,'bed');
% % bm_s =  bedmap2_interp(Xi,Yi,'surface');

bm_b =  bedmachine_interp('bed',Xi,Yi);
bm_s =  bedmachine_interp('surface',Xi,Yi);

if(runType == 3)
    %Golledge runs
    load Golledge21_GRL_T2_thick_22mar23_v2_Paul.mat
    [xgrid,ygrid] = meshgrid(is2xvec,is2yvec);
    s_interp = griddedInterpolant(xgrid',ygrid',surf_interp(:,:,thin_m)','linear','nearest');
    b_interp = griddedInterpolant(xgrid',ygrid',bed_interp(:,:,thin_m)','linear','nearest');
    goll_b =  b_interp(Xi,Yi);
    goll_s =  s_interp(Xi,Yi);

    clear is2xvec is2yvec xgrid ygrid;

    smoothbed = goll_b;%imgaussfilt(bm_b,2e3*2/dx);
    smoothsurf = imgaussfilt(goll_s,10e3/dx);
else
    smoothbed = bm_b;%imgaussfilt(bm_b,2e3*2/dx);
    smoothsurf = imgaussfilt(bm_s,10e3/dx);
end

% Smoothing
% Numerator is the window we're smoothing over in [m], spacing of these grids
% is actually dx/2 for bm_X grids hence the extra "*2".
% 

smooth_bm_bed = bm_b;%imgaussfilt(bm_b,2e3*2/dx);
smooth_bm_surf = imgaussfilt(bm_s,10e3/dx);

% smoothbed = sgolayfilt(bm_b,2,2*floor(10e3/dx)+1);
% smoothsurf = sgolayfilt(bm_s,2,2*floor(10e3/dx)+1);

% rectify rock above ice issue, force that ice is non-zero thickness everywhere
smoothsurf(smoothbed > smoothsurf) = smoothbed(smoothbed > smoothsurf) + 1; %Pe and Br = 0 result in NaNs


%% Build bed and surf, correct for thinning and floatation
if(runType == 2)
    load ATL15_dhdt.mat
    [xgrid,ygrid] = meshgrid(xvec,yvec);
    dhdt_interp = griddedInterpolant(xgrid',ygrid',dhdt','linear','nearest');
    clear xvec yvec xgrid ygrid;
end

h_real =@(x,y) interp2(xi,yi,bm_s-bm_b,x,y);
rock_mask =@(x,y) interp2(xi,yi,rock,x,y,'nearest');
h_bm_b =@(x,y) interp2(xi,yi,smooth_bm_bed,x,y);
h_bm_s =@(x,y) interp2(xi,yi,smooth_bm_surf,x,y);
h_b_init =@(x,y) interp2(xi,yi,smoothbed,x,y);

if(runType == 1)
    h_s_init =@(x,y) interp2(xi,yi,smoothsurf,x,y) + thin_m;
elseif(runType == 2)
    h_s_init =@(x,y) interp2(xi,yi,smoothsurf,x,y) - thin_m.* dhdt_interp(x,y); % minus to go back in time
elseif(runType == 3)
    h_s_init =@(x,y) interp2(xi,yi,smoothsurf,x,y); %Case where thin_m controls case of Golledge runs
else
    error("unknown runType")
end

phi_init =@(x,y) rho/rho_w*h_s_init(x,y) + (rho_w-rho)/rho_w*h_b_init(x,y); %hydropotential per unit water weight
clear bm_b bm_s;
h_init =@(x,y) subplus(h_s_init(x,y) - h_b_init(x,y)-1)+1; %h: smoothed ice thickness [m]
h_bm =@(x,y) subplus(h_bm_s(x,y) - h_bm_b(x,y)-1)+1;
h = h_init(xy(:,1),xy(:,2));

%% Define a few globals vars
phi_max = max(max(phi_init(xy(:,1),xy(:,2))));
phi_min = min(min(phi_init(xy(:,1),xy(:,2))));

%% Mask Between rock/sed (0 is rock, 1 is sed)

rockSed = double(Yi > -(Xi + 5e5)/1.7); %divide between rock/sed midstream
rockSed2 = heaviside(Yi + (Xi + 5e5)/1.7);
rockSedMask = griddedInterpolant(Xi',Yi',rockSed','nearest');   

%% Create vectors of bed/surface for numerical solving
h_s = h_s_init(xy(:,1),xy(:,2));
h_b = h_b_init(xy(:,1),xy(:,2));
