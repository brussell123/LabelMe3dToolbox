function LM2VRML(annotation,mesh,vrmlfile,vrmlfolder)
% Creates VRML 3D file of LabelMe3D scene.
%
% Inputs:
% annotation - LabelMe annotation structure.
% mesh - 3D triangular mesh structure.
% vrmlfile - VRML filename (e.g. 'example.wrl')
% vrmlfolder - Folder to dump textures to

% VRML Parameters:
fieldOfView = 0.78;

if nargin < 4
  vrmlfolder = pwd;
end

% Save textures 
relativeFolder = strrep(vrmlfile, '.wrl', '');
if isfield(mesh,'textures')
  mkdir(fullfile(vrmlfolder,relativeFolder));
  for ii = 1:length(mesh.textures)
    if ~isempty(mesh.textures{ii})
      texturefilename{ii} = sprintf('object_%03d.jpg', ii);
      imwrite(mesh.textures{ii},fullfile(vrmlfolder,relativeFolder,texturefilename{ii}));
    end
  end
else
  colors = hsv(length(annotation.object));
end

% Get initial camera position:
P = getCameraMatrix(annotation,'RH');
[K,R,initialCamera] = decomposeP(P);

% Dump VRML file:
fp = fopen(fullfile(vrmlfolder,vrmlfile),'w');
fprintf(fp,'#VRML V2.0 utf8\n\nNavigationInfo {\n  headlight TRUE\n  type [\"EXAMINE\", \"ANY\"]\n}\n\n');
fprintf(fp,'Viewpoint {\n  position %.2f %.2f %.2f\n  orientation 0 1 0 3.14159\n  fieldOfView %.2f\n  description "Original"\n}\n\n',initialCamera(1), initialCamera(2), initialCamera(3), fieldOfView);
fprintf(fp,'DEF Back1 Background {\n  groundColor [1 1 1]\n  skyColor [1 1 1]\n}\n\n');

for ii = 1:length(annotation.object)
  isValidObject = any(mesh.objndx==ii);
  if isValidObject
    fprintf(fp,'# %s\nShape{\n',annotation.object(ii).name);
    if isfield(mesh,'textures')
      fprintf(fp,'  appearance Appearance {\n    texture ImageTexture { url \"%s\" }\n  }\n', fullfile(relativeFolder,texturefilename{ii}));
    else
      fprintf(fp,'  appearance Appearance {\n    material Material {\n      diffuseColor %.2f %.2f %.2f\n    }\n  }\n',colors(ii,:));
    end

    % Get object mesh and texture coordinates:
    mesh_i = getObjectMesh(mesh,ii);

    fprintf(fp,'  geometry IndexedFaceSet {\n    coord Coordinate {\n      point [\n');
    for jj = 1:size(mesh_i.vertices,2)
      fprintf(fp,'        %.4f %.4f %.4f,\n',mesh_i.vertices(:,jj));
    end
    fprintf(fp,'      ]\n    }\n');
    
    %% Get coordIndex:
    fprintf(fp,'    coordIndex [\n');
    for j = 1:size(mesh_i.faces,2)
      fprintf(fp,'      %d %d %d -1,\n',mesh_i.faces(:,j)-1);
    end
    fprintf(fp,'    ]\n');
    
    if isfield(mesh,'textures')
      %% Get texture coordinates:
      fprintf(fp,'    texCoord TextureCoordinate {\n      point [\n');
      for jj = 1:size(mesh_i.tx,2)
        fprintf(fp,'        %.4f %.4f,\n',mesh_i.tx(1,jj)/size(mesh_i.textures{ii},2),1-mesh_i.tx(2,jj)/size(mesh_i.textures{ii},1));
      end
      fprintf(fp,'      ]\n');
      fprintf(fp,'    }\n');
      
      %% Get texture coordinate indices:
      fprintf(fp,'    texCoordIndex [\n');
      for j = 1:size(mesh_i.faces,2)
        fprintf(fp,'      %d %d %d -1,\n',mesh_i.faces(:,j)-1);
      end
      fprintf(fp,'    ]\n');
    end
    
    fprintf(fp,'    solid FALSE\n');
    fprintf(fp,'  }\n');
    fprintf(fp,'}\n\n');
  end
end
fclose(fp);
