function [E] = depthIntEnhancement(T_z,b_0,nn,dz)
%[A] = depthIntEnhancement(T_z, H, dz)
%returns depth averaged value of stiffness enhancement, following Minchew et al
%no lubrication of crystals is included, only temp viscous dependence.
    E = trapz(calcAfromT(T_z).^(-1/nn),2)/b_0*dz;
end