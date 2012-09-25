function annotation = InflateSceneNew(annotation)
% annotation = InflateScene(annotation)
%
% Given extracted image information, inflate scene.
%
% Inputs:
% annotation - LabelMe annotation structure.
% 
% Outputs:
% annotation - LabelMe annotation with "mesh" field for each object

% $$$ % Insert ground plane:
% $$$ annotation = insertGroundPlane(annotation);

% Get indices of non-deleted objects:
notDeleted = find(~isdeleted(annotation))';

Nobjects = length(annotation.object); % Number of objects

if isfield(annotation.object,'mesh')
  annotation.object = rmfield(annotation.object,'mesh');
end

for i = notDeleted%1:Nobjects
  polyType = annotation.object(i).polygon.polyType;
  if strcmp(polyType,'ground') || strcmp(polyType,'standing')
    [X,Y,Z,TX,TY,tri,boundary,CptsNdx,valid] = getObject3D(annotation,i,'trimesh');
    if valid
      annotation.object(i).mesh.v3D = [X; Y; Z];
      annotation.object(i).mesh.v2D = [TX; TY];
      annotation.object(i).mesh.tri = tri;
      annotation.object(i).mesh.boundary = boundary;
% $$$     annotation.object(i).mesh.type = polytype;
% $$$     if strcmp(polytype,'foldedcard')
      if strcmp(polyType,'standing')
        annotation.object(i).mesh.CptsNdx = CptsNdx;
      end
    end
  end
end
