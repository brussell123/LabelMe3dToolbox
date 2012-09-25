function [X,Y,Z,x,y,CptsNdx,tri,boundary] = standingFoldedCardWithMesh(x,y,CptsNdx,ncols,nrows,CAM_H,Hy,F)
% [X,Y,Z,x,y,CptsNdx,tri,boundary] = standingFoldedCardWithMesh(x,y,CptsNdx,ncols,nrows,CAM_H,Hy,F)
%
% This function estimates the depth of a shape that is like a folded card
% directly standing on the ground. For objects with a unique contact point,
% it assumes that the object is a flat surface parallel to the camera
% plane.
%
% Input:
%    x,y: image polygon coordinates
%    CptsNdx = [0 0 1 1 1 0 ...] binary indicator vector for contact
%                                points, same length than x and y
%    ncols, nrows: image size
%    CAM_H, Hy, F:  camera parameters
%
% Output: meshes
%    X,Y,Z: world coordinates of the polygon
%    x,y: new polygon coordinates (with inserted points)
%    CptsNdx: indicator of contact points
%    tri: triangular mesh indices
%    boundary: indices to points in the boundary

% Parameters:
npts = 11;

% Get contact points:
[x,y,CptsNdx] = addVerticalPoint(x,y,unique(x),CptsNdx);
j = find(CptsNdx);
xc = x(j);
yc = y(j);

% Need to have unique xc, else interp1 does not work.  Need to handle
% this in a better way.
[xc,i] = unique(xc);
yc = yc(i);

if nargout>6
  % Create mesh x,y:
  [x,y,tri,boundary] = getMesh(x,y,npts);
end

% Convert image coordinates into world coordinates:
[X,Y,Z] = Image2WorldCoords(x,y,xc,yc,[nrows ncols],CAM_H,Hy,F);
