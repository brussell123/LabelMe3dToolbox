% This script demonstrates the main functionalities of the LabelMe3D
% toolbox.  The functionalities depend on the LabelMe toolbox, which can
% be downloaded from here:
%
% http://labelme.csail.mit.edu/LabelMeToolbox/index.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 0: Compile mex files in the toolbox
compile;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Change this paths to point to the desired location of the
% LabelMe3D database
HOMEDATABASE = './LM3D_database';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Add required toolboxes to path:
addpath '~/work/MatlabLibraries/LabelMeToolbox';
addpath '~/work/MatlabLibraries/LabelMe3dToolbox';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define needed paths for the database:
HOMEIMAGES = fullfile(HOMEDATABASE,'Images');
HOMEANNOTATIONS = fullfile(HOMEDATABASE,'Annotations');
CACHE_DIR = fullfile(HOMEDATABASE,'cache'); % Location of any 3D outputs

% Download data:
LM3Dinstall(HOMEIMAGES,HOMEANNOTATIONS);

% Read data into Matlab structure.  The LabelMe annotation structure for
% the ith image is stored as "DB(i).annotation".
DB = LM3Ddatabase(HOMEANNOTATIONS);

% Read an annotation and image directly from the LabelMe server:
HOMEANNOTATIONS = 'http://labelme.csail.mit.edu/Annotations';
HOMEIMAGES = 'http://labelme.csail.mit.edu/Images';
folderlist = {'05june05_static_street_boston'};
filelist = {'p1010736.jpg'};
DB = LM3Ddatabase(HOMEANNOTATIONS,HOMEIMAGES,folderlist,filelist);

% Get an image and annotation:
i = 1;
annotation = DB(i).annotation;
img = LMimread(DB,i,HOMEIMAGES);

% Plot 3D annotations:
figure;
LMplot3Dscene(annotation);

% Create depth maps:

% Generate textured mesh:
mesh = LM3DgetTexturedMesh(annotation,img);

% Generate VRML file:
vrmlfile = strrep(annotation.filename,'.jpg','.wrl');
vrmlfolder = fullfile(CACHE_DIR,'vrml');
if ~exist(vrmlfolder,'dir')
  mkdir(vrmlfolder);
end
LM2VRML(annotation,mesh,vrmlfile,vrmlfolder);
display(sprintf('Produced a VRML file here: %s',fullfile(vrmlfolder,vrmlfile)));

% Sanity check that 3D points project to 2D points:
err = [];
for i = 1:length(annotation.object)
  if isfield(annotation.object(i).world3d,'polygon3d')
    X = getLMpolygon3D(annotation.object(i).world3d.polygon3d); % Get 3D points
    [x,y] = Project3D2D(X,annotation); % Project 3D points to 2D
    [xx,yy] = getLMpolygon(annotation.object(i).polygon); % Original 2D polygon
    err = [err [abs(x-xx'); abs(y-yy')]]; % Record error
  end
end
display(sprintf('Mean error: %f',mean(err(:))));

% Plot labeled objects and horizon line:
[xh,yh] = getHorizonLine(annotation);
figure;
LMplot(annotation,img);
hold on;
plot(xh,yh,'r');

% Run LabelMe3D algorithm.
% WARNING: this will over-write your XML annotations!
HOMEANNOTATIONS = fullfile(HOMEDATABASE,'Annotations');
LM3Dgenerate3D(HOMEANNOTATIONS,HOMEIMAGES);
