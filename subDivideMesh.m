function mesh = subDivideMesh(mesh,levels)

if nargin < 2
  levels = 1;
end

for i = 1:levels
  Nvertices = size(mesh.vertices,2);
  Nfaces = size(mesh.faces,2);
  v1 = [(mesh.vertices(1,mesh.faces(1,:))+mesh.vertices(1,mesh.faces(2,:)))/2; ...
        (mesh.vertices(2,mesh.faces(1,:))+mesh.vertices(2,mesh.faces(2,:)))/2; ...
        (mesh.vertices(3,mesh.faces(1,:))+mesh.vertices(3,mesh.faces(2,:)))/2];
  
  v2 = [(mesh.vertices(1,mesh.faces(2,:))+mesh.vertices(1,mesh.faces(3,:)))/2; ...
        (mesh.vertices(2,mesh.faces(2,:))+mesh.vertices(2,mesh.faces(3,:)))/2; ...
        (mesh.vertices(3,mesh.faces(2,:))+mesh.vertices(3,mesh.faces(3,:)))/2];
  
  v3 = [(mesh.vertices(1,mesh.faces(3,:))+mesh.vertices(1,mesh.faces(1,:)))/2; ...
        (mesh.vertices(2,mesh.faces(3,:))+mesh.vertices(2,mesh.faces(1,:)))/2; ...
        (mesh.vertices(3,mesh.faces(3,:))+mesh.vertices(3,mesh.faces(1,:)))/2];
  
  [v,ii,jj] = unique([v1 v2 v3]','rows' );
  
  m1 = int32(jj(1:Nfaces)+Nvertices)';
  m2 = int32(jj(Nfaces+1:2*Nfaces)+Nvertices)';
  m3 = int32(jj(2*Nfaces+1:3*Nfaces)+Nvertices)';
  
  tri1 = [mesh.faces(1,:); m1; m3];
  tri2 = [mesh.faces(2,:); m2; m1];
  tri3 = [m1; m2; m3];
  tri4 = [m2; mesh.faces(3,:); m3];
  
  mesh.vertices = [mesh.vertices v']; % the new vertices
  mesh.faces = [tri1 tri2 tri3 tri4]; % the new faces
  mesh.objndx = repmat(mesh.objndx,1,4);
end

% $$$ [vertices,faces,ndx] = subDivideMesh_helper(mesh.vertices,mesh.faces);
% $$$ 
% $$$ mesh.vertices = vertices;
% $$$ mesh.faces = faces;
% $$$ mesh.objndx = mesh.objndx(ndx);

% Invalidate texture coordinates:
if isfield(mesh,'tx')
  mesh = rmfield(mesh,'tx');
end


return;

i = 6;
mesh_i = getObjectMesh(mesh,i);


addpath '~/work/MatlabLibraries/LabelMeToolbox';
addpath '~/work/MatlabLibraries/LabelMe3dToolbox';
addpath './NewCodeToInclude';

load foo.mat;

mesh_out = subDivideMesh(mesh_i);

figure;
trisurf(mesh_out.faces',mesh_out.vertices(1,:),mesh_out.vertices(2,:),mesh_out.vertices(3,:));
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
cameratoolbar('setmodeGUI','orbit');

figure;
trisurf(mesh_i.faces',mesh_i.vertices(1,:),mesh_i.vertices(2,:),mesh_i.vertices(3,:));
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
cameratoolbar('setmodeGUI','orbit');
