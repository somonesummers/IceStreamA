function [out] = cleanNan(x,y,in)
% function [out] = cleanNan(x,y,in)
% Interpolates NAN values, works for all size arrays
    if (numel(x) ~= numel(y) || numel(in) ~= numel(x))
        error('inconsistent number of elements in arrays')
    end
    inLin = in(:);
    xLin = x(:);
    yLin = y(:);
    inFill = scatteredInterpolant(xLin(~isnan(inLin)),yLin(~isnan(inLin)),inLin(~isnan(inLin)));
    inLin(isnan(inLin)) = inFill(xLin(isnan(inLin)),yLin(isnan(inLin)));
    out = reshape(inLin, size(in));
end