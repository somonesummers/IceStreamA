function [tau_c] = defineTau(str,x0)
% [tau_c] = defineTau(str,[x0])
% returns tau and string for requested string-name of tau scenario. 
    opt = false;
    if(nargin == 2)
        opt = true;
    end
    if(str == "ISSM")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.38;%1.38;
        floor = 0e3;
    end
    warning('off','MATLAB:imagesci:netcdf:fillValueTypeMismatch'); % try not to warn here
    if(ismac)
        xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
        yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
        tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    else
        xi   = ncread("/home/groups/jsuckale/psummers/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
        yi   = ncread("/home/groups/jsuckale/psummers/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
        tau  = ncread("/home/groups/jsuckale/psummers/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    end
    warning('on','MATLAB:imagesci:netcdf:fillValueTypeMismatch'); % don't warn here
%     [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant({xi - 3072000,yi - 3072000},tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-floor)+floor)*scale;
    elseif(str == "ISSM Tuned")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.05; %1.2
        floor = 0e3;
    end
    load tau_tuned.mat;
    load gridSiple10000.mat;
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_tuned,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-floor)+floor)*scale;
    elseif(str == "ISSM Shift")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        stag = x0(1);
        moving = x0(2);
    else
        stag = 1.3;%1.55
        moving = 1.0;%.7
    end
    load tauShiftable.mat;
    load gridSiple1000.mat;
    
    tau_shift = zeros(size(tau_ISSM));
    tau_shift(spd<50) = tau_ISSM(spd<50).* stag; 
    tau_shift(spd>50) = tau_ISSM(spd>50).* moving; 
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_shift,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)));
    elseif(str == "ISSM Overburden")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        stag = x0(1);
        moving = x0(2);
    else
        stag = 1;%1.55
        moving = 1;%.7
    end
    load tauShiftable.mat;
    load gridSiple1000.mat;
    
    tau_shift = zeros(size(tau_ISSM));
    tau_shift(spd<50) = tau_ISSM(spd<50).* stag; 
    tau_shift(spd>50) = tau_ISSM(spd>50).* moving; 
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_shift,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)));
    elseif(str == "ISSM Shift3")%2.1 stag, .7 moving from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.00; %1.2
        floor = 0e3;
    end
    load tau_shift3.mat;
    load gridSiple1000.mat;
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_shift,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
       scale.* (subplus(uB(x,y)));
   elseif(str == "ISSM Shift4")  %1.7 stag, .7 moving from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.00; %1.2
        floor = 0e3;
    end
    load tau_shift4.mat;
    load gridSiple1000.mat;
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_shift,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
       scale.* (subplus(uB(x,y)));
   elseif(str == "ISSM Shift5")  %1.5 stag, .7 moving from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.00; %1.2
        floor = 0e3;
    end
    load tau_shift5.mat;
    load gridSiple1000.mat;
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_shift,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
       scale.* (subplus(uB(x,y)));
   elseif(str == "ISSM Shift6")  %1.3 stag, .7 moving from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.00; %1.2
        floor = 0e3;
    end
    load tau_shift6.mat;
    load gridSiple1000.mat;
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_shift,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
       scale.* (subplus(uB(x,y)));
    elseif(str == "PISM1")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 1.1;
        floor = 0;
    end
    xi   = ncread("~/Documents/MATLAB/ISSM/ARC_PISM1_ctrl/strbasemag_AIS_ARC_PISM1_ctrl.nc","x");
    yi   = ncread("~/Documents/MATLAB/ISSM/ARC_PISM1_ctrl/strbasemag_AIS_ARC_PISM1_ctrl.nc","y");
    tau  = ncread("~/Documents/MATLAB/ISSM/ARC_PISM1_ctrl/strbasemag_AIS_ARC_PISM1_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi,yi);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-floor)+floor)*scale;
    elseif(str == "PISM2")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 2;
    end
    xi   = ncread("~/Documents/MATLAB/ISSM/ARC_PISM2_ctrl/strbasemag_AIS_ARC_PISM2_ctrl.nc","x");
    yi   = ncread("~/Documents/MATLAB/ISSM/ARC_PISM2_ctrl/strbasemag_AIS_ARC_PISM2_ctrl.nc","y");
    tau  = ncread("~/Documents/MATLAB/ISSM/ARC_PISM2_ctrl/strbasemag_AIS_ARC_PISM2_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi,yi);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        subplus(uB(x,y))*scale;
    elseif(str == "ELMER")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 1.8;
    end
    xi   = ncread("~/Documents/MATLAB/ISSM/IGE_ELMER_ctrl/strbasemag_AIS_IGE_ELMER_ctrl.nc","x");
    yi   = ncread("~/Documents/MATLAB/ISSM/IGE_ELMER_ctrl/strbasemag_AIS_IGE_ELMER_ctrl.nc","y");
    tau  = ncread("~/Documents/MATLAB/ISSM/IGE_ELMER_ctrl/strbasemag_AIS_IGE_ELMER_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi,yi);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-10e3)+10e3)*scale;
    
    else
        error("Invalid Tau_c Scenario String")
    end
end
