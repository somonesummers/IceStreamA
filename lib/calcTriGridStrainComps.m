function [int1,int2] = calcTriGridStrainComps(u,v,xy,dx)
% tic
    % [int1, int2] = calcTrigridStrainComps(u,v,xy,dx) 
    % Calculates strain on triangular nodes, returns on intpolation
    % functions for strain1 and strain2 priciple strains
    
    % elements
    xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
    yi = (min(xy(:,2))-dx:dx:max(xy(:,2))+dx)';
    [Xi,Yi] = ndgrid(xi,yi);
    us = scatteredInterpolant(xy(:,1),xy(:,2),u);
    vs = scatteredInterpolant(xy(:,1),xy(:,2),v);

    %% grid to find strain rates
    ug = us(Xi,Yi); % [m/s]
    vg = vs(Xi,Yi); % [m/s]

    %% Smooth
%     ug = imgaussfilt(ug,1);
%     vg = imgaussfilt(vg,1);
    
    %% Calc Strain on grid
    [n,m] = size(ug);
    strain = zeros(n-1,m-1,2,2);
    strain(:,:,1,1) = diff(vg(1:end-1,:),1,2)/dx;
    strain(:,:,2,2) = diff(ug(:,1:end-1))/dx;
    strain(:,:,1,2) = .5 .* (diff(ug(1:end-1,:),1,2) + diff(vg(:,1:end-1)))/dx;
    strain(:,:,2,1) =  strain(:,:,1,2);
    
%%   Calc principle Shear Strains
   strain1 = zeros(n-1,m-1);
   strain2 = zeros(n-1,m-1);
     for i = 1:n-1   
         for j = 1:m-1
            M = squeeze(strain(i,j,:,:));
            if(sum(sum(isnan(M)))==0)
                [~,D] = eig(M);
                e = diag(D);
                strain1(i,j) = e(1);
                strain2(i,j) = e(2);
            else
                strain1(i,j) = 1e-21; %when nan values toss in near 0 value
                strain2(i,j) = 1e-21;
            end
         end
     end
     %% Interpolate back to triangular grid elements
     int1 = griddedInterpolant(Xi(2:end,2:end)-dx/2, Yi(2:end,2:end)-dx/2, strain1,'linear','nearest'); %Evaluate at middle of cell
     int2 = griddedInterpolant(Xi(2:end,2:end)-dx/2, Yi(2:end,2:end)-dx/2, strain2,'linear','nearest'); %Evaluate at middle of cell
%  disp("calcTriGridStrainComps: " + toc);
end