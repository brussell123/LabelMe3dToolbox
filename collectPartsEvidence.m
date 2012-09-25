function [evidence,objectClasses] = collectPartsEvidence(DB,objectClasses)
% [evidence,objectClasses] = collectPartsEvidence(DB)
% or
% evidence = collectPartsEvidence(annotation,objectClasses,validObjects)
%
% Collect evidence for attachment (part-of) relationships.
%
% Inputs:
% DB (annotation) - LabelMe DB structure.
% objectClasses - Cell array of object classes.
%
% Outputs:
% evidence(i).relativeOverlap - Relative overlap.
% evidence(i).relativeArea - Relative area.
% evidence(i).objNdx - Object class index.

if ~isfield(DB,'annotation')
  DB.annotation = DB;
end
if nargin < 2
  objectClasses = unique(lower(LMobjectnames(DB,'name')));
end

Nobjects = length(objectClasses);
Nimages = length(DB);
for i = 1:Nimages
  if isfield(DB(i).annotation,'object')
    if isfield(DB(i).annotation,'imagesize')
      nrows = DB(i).annotation.imagesize.nrows;
      ncols = DB(i).annotation.imagesize.ncols;
    else
      [ncols,nrows] = getaproximagesize(DB(i).annotation);
    end
    
    % Get object indices and relative area:
    Nobj = length(DB(i).annotation.object);
    objNdx = zeros(1,Nobj,'uint16');
    relativeArea = zeros(1,Nobj,'single');
    clear X Y;
    for j = 1:Nobj
      k = strmatch(strtrim(lower(DB(i).annotation.object(j).name)),objectClasses,'exact');
      if ~isempty(k)
        objNdx(j) = k;
      end
      [X{j},Y{j}] = getLMpolygon(DB(i).annotation.object(j).polygon);
      relativeArea(j) = min(max(polyarea(X{j},Y{j})/ncols/nrows,0),1);
    end

    % Collect relative overlaps:
    R = single(eye(Nobj));
    for m = 1:Nobj-1
      for n = m+1:Nobj
        [Rm,Rn] = relativeOverlap(X{m},Y{m},X{n},Y{n});
        R(m,n) = Rm;
        R(n,m) = Rn;
      end
    end
    
    % Store evidence into structure:
    evidence(i).relativeOverlap = R;
    evidence(i).relativeArea = relativeArea;
    evidence(i).objNdx = objNdx;
  else
    evidence(i).relativeOverlap = [];
    evidence(i).relativeArea = [];
    evidence(i).objNdx = [];
  end
end
