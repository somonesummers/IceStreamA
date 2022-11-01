% Used to initialize the input files for main runner

%% Color Packages
load Dawn.mat;
icey = cbrewer('div','BrBG',48);

%% Load Grid
load("grids/"+mapFile);

%% Load Lake Locations, make polygons in polar stereographic projection
% Sl = shaperead('tc-11-451-2017-supplement/Thw_lakes_supplemental_data/Thw_lakes_outlines.shp');
% Sn = size(Sl,1);
% for i = 1:Sn
%     [x,y] = ll2ps(Sl(i).Y(~isnan(Sl(i).Y)),Sl(i).X(~isnan(Sl(i).X)));
%     Sl(i).X = x;
%     Sl(i).Y = y;
%     clear x y
% end
% clear i

%% Get speed from measures for BCs
[spd_BC_u_se, spd_BC_v_se] = measures_interp('velocity',xy(se_bound,1),xy(se_bound,2));
[spd_BC_u_ne, spd_BC_v_ne] = measures_interp('velocity',xy(ne_bound,1),xy(ne_bound,2));
[spd_BC_u_sw, spd_BC_v_sw] = measures_interp('velocity',xy(sw_bound,1),xy(sw_bound,2));
[spd_BC_u_nw, spd_BC_v_nw] = measures_interp('velocity',xy(nw_bound,1),xy(nw_bound,2));

if(sum(isnan(spd_BC_v_se)) > 0)
    warning("NaN values in SE boundary condition: " + sum(isnan(spd_BC_v_se))); 
    spd_BC_u_se(isnan(spd_BC_u_se)) = 0;
    spd_BC_v_se(isnan(spd_BC_v_se)) = 0;
end
if(sum(isnan(spd_BC_v_ne)) > 0)
    warning("NaN values in NE boundary condition: " + sum(isnan(spd_BC_v_ne))); 
    spd_BC_u_ne(isnan(spd_BC_u_ne)) = 0;
    spd_BC_v_ne(isnan(spd_BC_v_ne)) = 0;
end
if(sum(isnan(spd_BC_v_sw)) > 0)
    warning("NaN values in SW boundary condition: " + sum(isnan(spd_BC_v_sw)));
    spd_BC_u_sw(isnan(spd_BC_u_sw)) = 0;
    spd_BC_v_sw(isnan(spd_BC_v_sw)) = 0;	
end
if(sum(isnan(spd_BC_v_nw)) > 0)
    warning("NaN values in NW boundary condition: " + sum(isnan(spd_BC_v_nw)));
    spd_BC_u_nw(isnan(spd_BC_u_nw)) = 0;
    spd_BC_v_nw(isnan(spd_BC_v_nw)) = 0;
end

%% Load accumlation and surface temp data
[Acc, T_s] = loadALBMAP();
T = T_s(xy(:,1),xy(:,2));
