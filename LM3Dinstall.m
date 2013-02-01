function LM3Dinstall(destImages, destAnnotations, folder)
%
% LM3Dinstall(destImages, destAnnotations)
%
% To download a single folder add a third argument:
% LM3Dinstall(destImages, destAnnotations, folder)

HOMEANNOTATIONS = 'http://labelme.csail.mit.edu/Annotations3D';
HOMEIMAGES = 'http://labelme.csail.mit.edu/Images';

% Build index
disp('Reading folder list')
if nargin == 3
    folder = folderlist(HOMEANNOTATIONS, folder);
else
    folder = folderlist(HOMEANNOTATIONS);
end

% create list of folders to copy, copy and rename folders in XML:
LMinstall(folder, destImages, destAnnotations, HOMEIMAGES, HOMEANNOTATIONS)

% make sure folder names are coherent
ADMINrenamefolder(destAnnotations);



