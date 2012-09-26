function caught_errs = LM3Dgenerate3D(HOMEANNOTATIONS,HOMEIMAGES)
% Runs LabelMe3D algorithm over annotations folder to generate 3D
% outputs.
%
% WARNING: This function will overwrite your XML annotations.
%
% Example: Run over all XML annotations inside HOMEANNOTATIONS folder:
%
% HOMEANNOTATIONS = '/path/to/Annotations';
% HOMEIMAGES = '/path/to/Images';
% LM3Dgenerate3D(HOMEANNOTATIONS,HOMEIMAGES)

% LabelMe3D structures:
load Params3D.mat;

% Make sure that folder names are consistent:
ADMINrenamefolder(HOMEANNOTATIONS);

% Get list of folders:
folders = strread(genpath(HOMEANNOTATIONS),'%s','delimiter',pathsep);

caught_errs = [];
for i = 1:length(folders)
  % Get list of XML files in folder:
  files = dir(fullfile(folders{i},'*.xml'));
  for j = 1:length(files)
    fname = fullfile(folders{i},files(j).name);
    try
      [annotation,img] = LMread(fname,HOMEIMAGES);
      imageSize = [size(img,1) size(img,2)];
      [annotation,valid] = Recover3DSceneComponents(annotation,Params3D,imageSize);
      if valid
        % Write XML file:
        writeXML3D(annotation,fname);
      end
    catch me
      display(sprintf('%s failed; continuing...',fname));
      if isempty(caught_errs)
        caught_errs.fname = fname;
      else
        caught_errs(end+1).fname = fname;
      end
      caught_errs(end).me = me;
    end
  end
end

return;


addpath '~/work/MatlabLibraries/LabelMeToolbox';
addpath '~/work/MatlabLibraries/LabelMe3dToolbox';

HOMEIMAGES = '/Users/brussell/work/Datasets/LabelMe3D/Images';
HOMEANNOTATIONS = '/Users/brussell/work/Datasets/LabelMe3D/Annotations';

caught_errs = LM3Dgenerate3D(HOMEANNOTATIONS,HOMEIMAGES);
