function mesh = getObjectMesh(mesh,i)
% Retrieve mesh for object.
%
% Inputs:
% mesh - Mesh structure.
% i - Object index.
%
% Outputs:
% mesh

% Get indices of faces belonging to object:
nFaces = find(mesh.objndx==i);

% Get indices of vertices belonging to object:
nVertices = unique(mesh.faces(:,nFaces));

% Get mapping from old vertices to new vertices:
mapVertices = zeros(1,size(mesh.vertices,2),'int32');
mapVertices(nVertices) = [1:length(nVertices)];

% Set output:
mesh.vertices = mesh.vertices(:,nVertices);
mesh.faces = mapVertices(mesh.faces(:,nFaces));
mesh.objndx = mesh.objndx(nFaces);
if isfield(mesh,'textures')
  texture = mesh.textures{i};
  mesh.textures = cell(1,length(mesh.textures));
  mesh.textures{i} = texture;
  mesh.tx = mesh.tx(:,nVertices);
end
