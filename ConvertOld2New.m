function annotation = ConvertOld2New(annotation,validObjects)
% Convert XML file to new format (this needs to be removed at some point)

for i = validObjects%1:length(annotation.object)
  annotation.object(i).world3d.stale = '0';
  annotation.object(i).world3d.type = '';

  % annotation.object.polygon.polyType -> annotation.object.world3d.type
  if isfield(annotation.object(i),'polygon') && isfield(annotation.object(i).polygon,'polyType')
    switch annotation.object(i).polygon.polyType
     case 'standing'
      annotation.object(i).world3d.type = 'standingplanes';
     case 'part'
      annotation.object(i).world3d.type = 'part';
     case 'ground'
      annotation.object(i).world3d.type = 'groundplane';
     otherwise
      annotation.object(i).world3d.type = 'none';
    end
    annotation.object(i).polygon = rmfield(annotation.object(i).polygon,'polyType');
  end

  % annotation.object.polygon.contact -> annotation.object.world3d.contact
  if isfield(annotation.object(i),'polygon') && isfield(annotation.object(i).polygon,'contact')
    annotation.object(i).world3d.contact = annotation.object(i).polygon.contact;
    annotation.object(i).polygon = rmfield(annotation.object(i).polygon,'contact');
  end

  % annotation.object.partof -> annotation.object.world3d.partof
  if isfield(annotation.object(i),'partof') && ~isempty(annotation.object(i).partof)
    annotation.object(i).world3d.parentid = annotation.object(i).partof;
  end
end
if isfield(annotation.object,'partof')
  annotation.object = rmfield(annotation.object,'partof');
end

% Convert image size to strings:
annotation.imagesize.nrows = num2str(annotation.imagesize.nrows);
annotation.imagesize.ncols = num2str(annotation.imagesize.ncols);
