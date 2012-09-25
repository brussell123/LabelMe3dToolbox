function [X,Y,Z,x,y,tri,boundary] = groundObjectWithMesh(x,y,ncols,nrows,CAM_H,Hy,F)
% [X,Y,Z,TX,TY,tri,boundary] = groundObjectWithMesh(x,y,ncols,nrows,CAM_H,Hy,F)
%
% This function estimates the depth of an horizontal shape lying on the ground 
% directly standing on the ground.
% Input:
%    x,y: image polygon coordinates
%    ncols, nrows: image size
%    CAM_H, Hy, F:  camera parameters
%
% Output:
%    X,Y,Z: world coordinates of the polygon
%    TX, TY: image coordinates, for the texture mesh
%    tri: triangular mesh indices
%    boundary: indices to points in the boundary

% Parameters:
npts = 11;

if nargout>5
  % Create mesh x,y:
  [x,y,tri,boundary] = getMesh(x,y,npts);
end

% Convert image coordinates into world coordinates:
[X,Y,Z] = Image2WorldCoords(x,y,[],[],[nrows ncols],CAM_H,Hy,F);
