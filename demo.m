% This script demonstrates the main functionalities of the LabelMe3D
% toolbox.  The functionalities depend on the LabelMe toolbox, which can
% be downloaded from here:
%
% http://labelme.csail.mit.edu/LabelMeToolbox/index.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: Compile mex files in the LabelMe3D toolbox
compile;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: Add required toolboxes to the Matlab path (change the following
% lines to point to the locations of the LabelMe and LabelMe3D toolboxes)
addpath '~/work/MatlabLibraries/LabelMeToolbox'; % LabelMe toolbox
addpath '~/work/MatlabLibraries/LabelMe3dToolbox'; % LabelMe3D toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOWNLOAD DATA, OPTION 1: Download the entire database onto your machine.

% Change this path to point to the desired location of the LabelMe3D database:
HOMEDATABASE = './LabelMe3D_database';

% Define needed paths for the database:
HOMEANNOTATIONS3D = fullfile(HOMEDATABASE,'Annotations3D');
HOMEIMAGES = fullfile(HOMEDATABASE,'Images');

% Download data:
LM3Dinstall(HOMEIMAGES,HOMEANNOTATIONS3D);

% Read data into Matlab structure.  The LabelMe annotation structure for
% the ith image is stored as "DB(i).annotation".
DB = LM3Ddatabase(HOMEANNOTATIONS3D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOWNLOAD DATA, OPTION 2: Download only needed images.

% Read an annotation and image directly from the LabelMe server:
HOMEANNOTATIONS3D = 'http://labelme.csail.mit.edu/Annotations3D';
HOMEIMAGES = 'http://labelme.csail.mit.edu/Images';
folderlist = {'05june05_static_street_boston'};
filelist = {'p1010736.jpg'};
DB = LM3Ddatabase(HOMEANNOTATIONS3D,HOMEIMAGES,folderlist,filelist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA VISUALIZATION.  The following illustrates how to visualize the data.

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
vrmlfolder = fullfile('cache','vrml');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUNNING LABELME3D ALGORITHM, OPTION 1: Run over all downloaded annotations.

% WARNING: this will over-write your downloaded XML annotations!
HOMEANNOTATIONS3D = fullfile(HOMEDATABASE,'Annotations3D');
LM3Dgenerate3D(HOMEANNOTATIONS3D,HOMEIMAGES);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUNNING LABELME3D ALGORITHM, OPTION 2: Run over a single LabelMe annotation.

% Get LabelMe annotation:
HOMEANNOTATIONS = 'http://labelme.csail.mit.edu/Annotations';
HOMEIMAGES = 'http://labelme.csail.mit.edu/Images';
folderlist = {'05june05_static_street_boston'};
filelist = {'p1010736.jpg'};
DB = LMdatabase(HOMEANNOTATIONS,HOMEIMAGES,folderlist,filelist);

% Run algorithm:
[annotation,valid] = Recover3DSceneComponents(DB.annotation);

if valid
  % Plot output:
  figure;
  LMplot3Dscene(annotation);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
