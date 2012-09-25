function [Xmax,Ymax,Zmax,valid] = getObjectDimensions(annotation,i)
% Inputs:
% annotation
% i
%
% Outputs:
% Xmax
% Ymax
% Zmax
% valid

Xmax = 0; Ymax = 0; Zmax = 0; valid = 0;
[X,Y,Z,mx,my,tri,boundary,CptsNdx,valid_i] = getObject3D(annotation,i,'poly');
if valid_i && ~isempty(X)
  Xmax = max(X)-min(X);
  Ymax = max(Y)-min(Y);
  Zmax = max(Z)-min(Z);
  valid = 1;
end
valid = logical(valid);

return;

notDeleted = find(~isdeleted(annotation))';
Nobjects = length(annotation.object);
Xmax = zeros(1,Nobjects);
Ymax = zeros(1,Nobjects);
Zmax = zeros(1,Nobjects);
valid = logical(zeros(1,Nobjects));
for i = notDeleted%1:Nobjects
  [X,Y,Z,mx,my,tri,boundary,CptsNdx,valid_i] = getObject3D(annotation,i,'poly');
  if valid_i && ~isempty(X)
    Xmax(i) = max(X)-min(X);
    Ymax(i) = max(Y)-min(Y);
    Zmax(i) = max(Z)-min(Z);
    valid(i) = 1;
  end
end

return;

Nobjects = length(annotation.object);
Xmax = zeros(1,Nobjects);
Ymax = zeros(1,Nobjects);
Zmax = zeros(1,Nobjects);
valid = logical(ones(1,Nobjects));
for i = 1:Nobjects
  j = getPartOfParents(annotation,i);
  j = j(end);
  if strcmp(annotation.object(i).polygon.polyType,'standing') || strcmp(annotation.object(i).polygon.polyType,'ground')
    [X,Y,Z] = getObject3D(annotation,i,'poly');
    Xmax(i) = max(X)-min(X);
    Ymax(i) = max(Y)-min(Y);
    Zmax(i) = max(Z)-min(Z);
  elseif strcmp(annotation.object(i).polygon.polyType,'part') && strcmp(annotation.object(j).polygon.polyType,'standing')
    [mx,my] = getLMpolygon(annotation.object(i).polygon);
    [X,Y,Z] = getObject3D(annotation,j,'mymesh',mx,my);
    Xmax(i) = max(X)-min(X);
    Ymax(i) = max(Y)-min(Y);
    Zmax(i) = max(Z)-min(Z);
  else
    valid(i) = 0;
  end
end
