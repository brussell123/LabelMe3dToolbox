function [Zmap,objectDepthMaps] = getDepthMap(annotation,img)
% Zmap = getDepthMap(annotation,img)
%
% Computes depth map given the computed 3D mesh.
%
% Inputs:
%   annotation - LabelMe annotation structure
%   img - Corresponding image
%
% Outputs:
%   Zmap - Depth map
%   objectDepthMaps - 3D matrix containing depths of the individual
%   objects.

MAX_DEPTH = 16777215;
Nobjects = length(annotation.object);
[nrows,ncols,dim] = size(img);

% Get camera parameters:
Hy = annotation.camera.Hy;
CAM_H = annotation.camera.CAM_H;
F = annotation.camera.F;

% Load predefined list of ground objects:
groundObjects = getListGroundObjects; 
groundObjects = sprintf('%s,', groundObjects{:}); groundObjects = groundObjects(1:end-1);
groundObjects = LMobjectindex(annotation, groundObjects, 'exact');

[mx,my] = meshgrid([1:ncols],[1:nrows]);
Zmap = MAX_DEPTH*ones(nrows,ncols,'uint32');
if nargout>=2
  objectDepthMaps = MAX_DEPTH*ones(nrows,ncols,Nobjects,'uint32');
end
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh)
    [X,Y] = getLMpolygon(annotation.object(i).polygon);
    v2D = annotation.object(i).mesh.v2D;
    v3D = annotation.object(i).mesh.v3D;
    n = find(inpolygon(mx,my,X,Y));
% $$$     z = uint32(griddata(v2D(1,:),v2D(2,:),v3D(3,:),mx(n),my(n)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This part should be abstracted since it also appears in
    % InflateScenesNew
    if ismember(i,groundObjects)
      [xx,yy,z] = Image2WorldCoords(mx(n),my(n),[],[],[nrows ncols],CAM_H,Hy,F);
      z = uint32(z);
    else
      % Get contact points:
      CptsNdx = [annotation.object(i).polygon.pt.contact]';
      [x,y,CptsNdx] = addVerticalPoint(X,Y,unique(X),CptsNdx);
      j = find(CptsNdx);
      xc = x(j);
      yc = y(j);
      
      % Need to have unique xc, else interp1 does not work.  Need to handle
      % this in a better way.
      [xc,j] = unique(xc);
      yc = yc(j);

      % Convert image coordinates into world coordinates:
      [xx,yy,z] = Image2WorldCoords(mx(n),my(n),xc,yc,[nrows ncols],CAM_H,Hy,F);
      z = uint32(z);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargout>=2
      objectDepthMaps(n+(i-1)*nrows*ncols) = z;
    end
    Zmap(n) = min(Zmap(n),z);
  end
end

if nargout==0
  plotDepthMap(Zmap);
end
