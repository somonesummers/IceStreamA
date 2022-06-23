function [tau_c] = defineTau(str,x0)
% [tau_c] = defineTau(str,h_s_init,h_b_init,phi_init,phi_max,phi_min,[x0])
% returns tau and string for requested string-name of tau scenario. 
    opt = false;
    if(nargin == 2)
        opt = true;
    end
    if(str == "Uniform") % Uniform Plastic Bed, for reference, not in paper
     if(opt)
        yield_base = x0;
    else
        yield_base = 100e3;
     end
    
    tau_c = @(x,y,u,v) norms([u,v],2,2)... 
       .*yield_base; 
    elseif(str == "VDV")  
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else 
        scale = 4;
        floor = 10e3;
    end
    a = load('bedDragDx5000smth10000Capped.mat');

    uB = griddedInterpolant(a.Xi',a.Yi',abs(a.bed'));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-floor)+floor)*scale;
    elseif(str == "ISSM")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 3;
        floor = 0e3;
    end
    xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
    yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
    tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
    [xx,yy] = ndgrid(xi - 3072000,yi - 3072000);
    uB = griddedInterpolant(xx,yy,tau(:,:,21));
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-floor)+floor)*scale;
    elseif(str == "ISSM Tuned")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
        floor = x0(2);
    else
        scale = 1.5; %1.2
        floor = 0e3;
    end
    load tau_tuned.mat;
    load gridSiple10000.mat;
    
    uB = scatteredInterpolant(xy(:,1),xy(:,2),tau_tuned,'natural');
    
    tau_c = @(x,y,u,v) norms([u,v],2,2) .*... %Plastic
        (subplus(uB(x,y)-floor)+floor)*scale;
    elseif(str == "PISM1")  % from https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.html
    if(opt)
        scale = x0(1);
    else 
        scale = 3;
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
