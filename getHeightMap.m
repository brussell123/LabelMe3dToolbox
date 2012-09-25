function [Ymap,objectHeightMaps] = getHeightMap(annotation,img,Zmap,objectDepthMaps)
% Zmap = getDepthMap(annotation,img)
%
% Computes depth map given the computed 3D mesh.
%
% Inputs:
%   annotation - LabelMe annotation structure
%   img - Corresponding image
%
% Outputs:
%   Ymap - Depth map
%   objectHeightMaps - 3D matrix containing depths of the individual
%   objects.

if nargin < 4
  [Zmap,objectDepthMaps] = getDepthMap(annotation,img);
end

MAX_HEIGHT = 16777215;
Nobjects = length(annotation.object);
[nrows,ncols,dim] = size(img);

[mx,my] = meshgrid([1:ncols],[1:nrows]);
Ymap = MAX_HEIGHT*ones(nrows,ncols,'uint32');
if nargout>=2
  objectHeightMaps = MAX_HEIGHT*ones(nrows,ncols,Nobjects,'uint32');
end
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh)
    [X,Y] = getLMpolygon(annotation.object(i).polygon);
    v2D = annotation.object(i).mesh.v2D;
    v3D = annotation.object(i).mesh.v3D;
    n = find(inpolygon(mx,my,X,Y));
    y = uint32(griddata(v2D(1,:),v2D(2,:),v3D(2,:),mx(n),my(n)));
    if nargout>=2
      objectHeightMaps(n+(i-1)*nrows*ncols) = y;
    end
    z = objectDepthMaps(:,:,i);
    j = find(Zmap(n)==z(n));
    Ymap(n(j)) = y(j);
  end
end

if nargout==0
  plotDepthMap(Ymap);
end
