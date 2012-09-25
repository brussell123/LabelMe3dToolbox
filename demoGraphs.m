HOMEANNOTATIONS = 'http://labelme.csail.mit.edu/Annotations';
HOMEIMAGES = 'http://labelme.csail.mit.edu/Images';

%addpath('.');
%addpath(genpath('./src'));

%% Create DB:
DB = LOADINDEX;

% remove poorly labeled images
labeled_area = LMlabeledarea(DB);
ImsNdxLabeled = find(labeled_area>=.9);
DB = DB(ImsNdxLabeled);
DB = LMaddtags(DB,'tags.txt');

%% Learn parts and support graphs:
objectClasses = LMobjectnames(DB,'name');

Ppart = partsGraph(DB, objectClasses);
Psup = supportGraph(DB, objectClasses, Ppart);
save graphsPartSuport objectClasses Ppart Psup DB

showPartsTree(Ppart, objectClasses, 'building', DB)
showPartsTree(Psup*3, objectClasses, 'road', DB)
