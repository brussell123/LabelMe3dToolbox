%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Change this paths to point to the desired location of the
% LabelMe3D database
HOMEDATABASE = './LM3D_database';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add required toolboxes to path:
addpath '~/work/MatlabLibraries/LabelMeToolbox';
addpath '~/work/MatlabLibraries/LabelMe3dToolbox';

% Define needed paths for the database:
HOMEIMAGES = fullfile(HOMEDATABASE,'Images');
HOMEANNOTATIONS = fullfile(HOMEDATABASE,'Annotations');
CACHE_DIR = fullfile(HOMEDATABASE,'cache'); % Location of any 3D outputs

% Download data:
LM3Dinstall(HOMEIMAGES,HOMEANNOTATIONS);

% Read data into Matlab structure:
DB = LMdatabase(HOMEANNOTATIONS);

% Get annotation:
i = 1;
annotation = DB(i).annotation;
img = LMimread(DB,i,HOMEIMAGES);

% Plot 3D annotations:

% Create depth maps:

% Generate textured mesh:
mesh = getSceneMesh(annotation);
mesh = subDivideMesh(mesh,2);
mesh = getTextures(mesh,annotation,img);

% Generate VRML file:
vrmlfile = strrep(annotation.filename,'.jpg','.wrl');
vrmlfolder = fullfile(CACHE_DIR,'vrml');
if ~exist(vrmlfolder,'dir')
  mkdir(vrmlfolder);
end
LM2VRML(annotation,mesh,vrmlfile,vrmlfolder);
display(sprintf('Produced a VRML file here: %s',fullfile(vrmlfolder,vrmlfile)));

% Run LabelMe3D algorithm.
% WARNING: this will over-write your XML annotations!
LM3Dgenerate3D(HOMEANNOTATIONS,HOMEIMAGES);
