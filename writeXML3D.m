function writeXML3D(annotation,fname)
% writeXML3D(annotation,fname)
%
% Write XML file with 3D information.
%
% Inputs:
% annotation - LabelMe annotation structure.
% fname - Filename.

if isfield(annotation,'object')
  if isfield(annotation.object,'mesh')
    annotation.object = rmfield(annotation.object,'mesh');
  end
  for ii = 1:length(annotation.object)
    if isfield(annotation.object(ii),'polygon')
      if isfield(annotation.object(ii).polygon,'edges')
        annotation.object(ii).polygon = rmfield(annotation.object(ii).polygon,'edges');
      end
      if isfield(annotation.object(ii).polygon,'pAttached')
        annotation.object(ii).polygon = rmfield(annotation.object(ii).polygon,'pAttached');
      end
      if isfield(annotation.object(ii).polygon,'pAttachedObject')
        annotation.object(ii).polygon = rmfield(annotation.object(ii).polygon,'pAttachedObject');
      end
    end
  end
end


v.annotation = annotation;
writeXML(fname,v);
