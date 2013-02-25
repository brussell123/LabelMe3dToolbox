function [annotation,valid] = Recover3DSceneComponents(annotation,Params3D,imageSize)
% Inputs:
% annotation - LabelMe annotation structure.
% Params3D - LabelMe3D parameters.
% imageSize - Image size [nrows ncols]
%
% Outputs:
% annotation - LabelMe annotation structure, augmented with 3D information.
% valid - Indicates whether output is valid.

valid = 0;
if (nargin < 3) && ~isfield(annotation,'imagesize')
  error('The annotation is missing the <imagesize> field.  Please pass in [nrows ncols] as third argument to Recover3DSceneComponents.m.');
  return;
end
if ~isfield(annotation,'imagesize')
  annotation.imagesize.nrows = imageSize(1);
  annotation.imagesize.ncols = imageSize(2);
end
if isstr(annotation.imagesize.nrows)
  annotation.imagesize.nrows = str2num(annotation.imagesize.nrows);
end
if isstr(annotation.imagesize.ncols)
  annotation.imagesize.ncols = str2num(annotation.imagesize.ncols);
end
imageSize = [annotation.imagesize.nrows annotation.imagesize.ncols];

if nargin < 2
  load Params3D.mat;
end

validObjects = [];
for i = 1:length(annotation.object)
  if ~isfield(annotation.object(i),'world3d') || ~isfield(annotation.object(i).world3d,'stale') || strcmp(annotation.object(i).world3d.stale,'1')
    validObjects(end+1) = i;
  end
end

% Add object IDs:
for i = 1:length(annotation.object)
  annotation.object(i).id = num2str(i-1);
end
maxID = annotation.object(end).id;

% Get valid object indices:
dd = isdeleted(annotation);
validObjects = intersect(validObjects,find(dd==0)');

display(sprintf('Computing geometry for %d objects.',length(validObjects)));

% Clean up object names:
annotation = myLMaddtags(annotation);

if length(validObjects)>0
  % Step 1: Infer polygon types (groundplane|standingplanes|part):
  annotation = inferPartsPosterior(annotation,Params3D.Ppart,Params3D.objectClasses,validObjects);
  
  % Step 2: Infer contact edges for "standingplanes" objects:
  annotation = inferEdgeTypesNew(annotation,validObjects);

  % Make sure all standing objects make contact with the ground:
  for i = validObjects
    if strcmp(annotation.object(i).polygon.polyType,'standing') && isempty(getLMcontact(annotation.object(i).polygon))
      annotation.object(i).polygon.polyType = '';
    end
  end
  
  % Convert XML file to new format (this needs to be removed at some point):
  annotation = ConvertOld2New(annotation,validObjects);
  
  % Step 3: Recover camera parameters:
% $$$   annotation = getviewpoint(annotation,Params3D.ObjectHeights);
  annotation = getviewpoint_ObjectHeights(annotation);
else
  display('Skipping part-of and edge type inference...');
end

% Get indices of non-deleted objects:
notDeleted = find(~isdeleted(annotation))';

for i = notDeleted
  if strcmp(annotation.object(i).world3d.type,'groundplane')
    valid = 1;
  end
end

% Restore original object names and remove "originalname" field:
for i = 1:length(annotation.object)
  annotation.object(i).name = annotation.object(i).originalname;
end
annotation.object = rmfield(annotation.object,'originalname');

% Remove points that have been added before:
for i = notDeleted
  switch annotation.object(i).world3d.type
   case {'standingplanes','part'}
    if isfield(annotation.object(i).world3d,'polygon3d')
      nAdded = logical(str2num(char({annotation.object(i).world3d.polygon3d.pt.added})));
      annotation.object(i).polygon.pt(nAdded) = [];
      annotation.object(i).world3d.polygon3d.pt(nAdded) = [];
    end
  end
end

% Inflate the 3D scene and store in annotation structure:
annotation = InflateSceneXML(annotation,notDeleted);
