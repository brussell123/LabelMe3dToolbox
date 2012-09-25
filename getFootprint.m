function ycc = getFootprint(xc,yc,x)
% ycc = getFootprint(xc,yc,x)
%
% Interpolate contact points.
%
% Inputs:
% xc - X coordinates of contact points
% yc - Y coordinates of contact points
% x - X coordinate to interpolate
%
% Output:
% ycc - Y coordinate of interpolated contact point

ycc = interp1(xc,yc,x,'linear','extrap');
