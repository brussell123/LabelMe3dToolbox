function [X,Y,Z,mx,my,tri,boundary,CptsNdx,valid] = getObject3D(annotation,i,meshtype,mx,my)
% Inputs:
% annotation - LabelMe annotation structure.
% i - Object index
% meshtype - 'grid' or 'trimesh' or 'poly' or 'mymesh'
%
% Outputs:
% X
% Y
% Z
% mx
% my
% tri
% boundary
% CptsNdx
% valid - Indicates whether 3D points are valid.

npts = 11; % Trimesh parameter

% Get camera matrix:
P = getCameraMatrix(annotation,'RH');
Hy = str2num(annotation.camera.Hy);

nrows = str2num(annotation.imagesize.nrows);
ncols = str2num(annotation.imagesize.ncols);

% Get polygon:
[x,y] = getLMpolygon(annotation.object(i).polygon);

% Get mesh for object:
switch meshtype
 case 'grid'
% $$$   [mx,my] = meshgrid([1:ncols],[1:nrows]);
% $$$   n = find(inpolygon(mx,my,x,y));
  bb = [max(1,floor(min(x))) min(ncols,ceil(max(x))) max(1,floor(min(y))) min(nrows,ceil(max(y)))];
  [mx,my] = meshgrid([bb(1):bb(2)],[bb(3):bb(4)]);
  n = inpolygon(mx,my,x,y);
  mx = mx(n);
  my = my(n);
  tri = [];
  boundary = [];
 case 'trimesh'
  [mx,my] = addVerticalPoint(x,y,unique(x));
  [mx,my,tri,boundary] = getMesh(mx,my,npts);
 case 'poly'
  [mx,my] = addVerticalPoint(x,y,unique(x));
  tri = [];
  boundary = [];
 case 'mymesh'
  tri = [];
  boundary = [];
 otherwise
  error('Invalid meshtype');
end

% Get 3D coordinates:
X = []; Y = []; Z = []; valid = 0;
CptsNdx = [];
switch annotation.object(i).polygon.polyType
 case 'ground'
  [X,Y,Z,valid] = Image2WorldCoords(mx,my,[],[],P);
 case 'standing'
  % Get contact points:
  [xc,yc] = getLMcontact(annotation.object(i).polygon);
  % Objects above the horizon line can not be inflated:
  if (max(y)>Hy) && ~isempty(xc)
    % Convert image coordinates into world coordinates:
    [X,Y,Z,valid] = Image2WorldCoords(mx,my,xc,yc,P);
  end
 case 'part'
  % Get parent object:
  parents = getPartOfParents(annotation,i);
  if ~isempty(parents)
    parentNdx = parents(end);
    switch annotation.object(parentNdx).polygon.polyType
     case 'ground'
      [X,Y,Z,valid] = Image2WorldCoords(mx,my,[],[],P);
     case 'standing'
      % Get contact points:
      [xc,yc] = getLMcontact(annotation.object(parentNdx).polygon);
      [x,y] = getLMpolygon(annotation.object(parentNdx).polygon);
      % Objects above the horizon line can not be inflated:
      if (max(y)>Hy) && ~isempty(xc)
        % Convert image coordinates into world coordinates:
        [X,Y,Z,valid] = Image2WorldCoords(mx,my,xc,yc,P);
      end
    end
  end
end
