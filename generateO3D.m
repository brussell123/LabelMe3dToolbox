function generateO3D(mesh,outfname)
% Inputs:
% mesh
% outfname

% Parameters:
scaling = 1/100; % Convert to meters

objNdx = unique(mesh.objndx);
Nobjects = length(objNdx);
fp = fopen(outfname,'w');
fprintf(fp,'%d',Nobjects);
for i = 1:Nobjects
  mesh_i = getObjectMesh(mesh,objNdx(i));
  mesh_i.tx(2,:) = -mesh_i.tx(2,:);
  fprintf(fp,',%d,%d,',length(mesh_i.vertices),length(mesh_i.faces));
  fprintf(fp,'%.4f,%.4f,%.4f,',mesh_i.vertices);
  fprintf(fp,'%.4f,%.4f,',mesh_i.tx);
  fprintf(fp,'%d,%d,%d,',mesh_i.faces);
end
fclose(fp);

return;

function generateO3D_old(annotation,img,outfname)
% Inputs:
% annotation
% img
% outfname

% Parameters:
scaling = 1/100; % Convert to meters

% Get meshes:
[X,Y,Z,TX,TY,tri,boundary,textures,objNames] = getTextures(annotation,img);

Nobjects = length(X);
fp = fopen(outfname,'w');
fprintf(fp,'%d',Nobjects);
for i = 1:Nobjects
  fprintf(fp,',%d,%d,',length(X{i}),length(tri{i}));
  for j = 1:length(X{i})
    fprintf(fp,'%.4f,%.4f,%.4f,',X{i}(j)*scaling,Y{i}(j)*scaling,-Z{i}(j)*scaling);
  end
  for j = 1:length(TX{i})
    fprintf(fp,'%.4f,%.4f,',TX{i}(j),1-TY{i}(j));
  end
  for j = 1:size(tri{i},1)-1
    fprintf(fp,'%d,%d,%d,',tri{i}(j,:));
  end
  j = size(tri{i},1);
  fprintf(fp,'%d,%d,%d',tri{i}(j,:));
end
fclose(fp);

return;

% Demo script:

addpath '~/work/MatlabLibraries/LabelMeToolbox';
addpath '~/work/MatlabLibraries/LabelMe3dToolbox';

fname_jpg = 'example.jpg';
fname_xml = 'example.out.xml';
fname_Z = 'example.Z.png';
img = imread(fname_jpg);
annotation = getfield(loadXML(fname_xml),'annotation');
Zmap = readZmap(fname_Z);
figure; LMplot(annotation,img); title('Image with annotations');
figure; plotZmap(Zmap); title('Depth map'); 

% Create object meshes (we need this for VRML and book outputs):
annotation = InflateSceneNew(annotation);

OUTDIR = './VRML';
% Create VRML file:
vrmlfilename  = 'example.wrl';
imwrite(img,fullfile(OUTDIR,'example.jpg'),'jpg','quality',95);
LM2VRMLfill(annotation,OUTDIR,vrmlfilename,'example.jpg',img);

o3dfilename = strrep(fname_jpg,'.jpg','_o3dmesh.txt');
generateO3D(annotation,img,o3dfilename);

% Parameters:
scaling = 1/100; % Convert to meters

% Get camera parameters:
fov = rad2deg(atan2(sqrt(size(img,1)^2+size(img,2)^2)/2,str2num(annotation.camera.F)));
cameraEye = [0 scaling*str2num(annotation.camera.CAM_H) 0];
sceneUp = [0 1 0];
C = 100; % Direction constant
cameraTarget = scaling*[0 str2num(annotation.camera.CAM_H)-C*(str2num(annotation.camera.Hy)-size(img,1)/2)/str2num(annotation.camera.F) -C];
