function annotation = inferPartsPosterior(annotation,Ppart,objectClasses,validObjects)
% annotation = inferPartsPosterior(annotation,Ppart,objectClasses,validObjects)
% Infer posterior of attachment relationship.
%
% Inputs:
% annotation - LabelMe annotation structure.
% Ppart - Learned part parameters.
% objectClasses - Cell array of object classes.
% validObjects - Vector of desired object indices to compute.
%
% Outputs:
% annotation.object(i).polygon.polyType - Polygon type ('part',
%                                         'standing','ground','')
% annotation.object(i).partof - ID of parent part-of polygon.
% annotation.object(i).polygon.pAttached - Probability that polygon is
%                                          attached.
% annotation.object(i).polygon.pAttachedObject - Probabilities for
%                                                part-of parents.

% Parameters:
threshPart = 0.4872; % Threshold for part-of relationship

if nargin < 4
  validObjects = 1:length(annotation.object);
end

notDeleted = (~isdeleted(annotation))';

% Collect parts evidence:
evidence = collectPartsEvidence(annotation,objectClasses);

% Get ground object indices:
groundObjects = getListGroundObjects; 
groundObjects = sprintf('%s,', groundObjects{:}); groundObjects = groundObjects(1:end-1);
groundObjects = LMobjectindex(annotation,groundObjects,'exact');

Nobjects = length(annotation.object);
for i = validObjects%1:Nobjects
  s = evidence.objNdx(i);

  if ismember(i,groundObjects)
    % This is a hack for now: ground object edges are automatically
    % considered attached:
    annotation.object(i).polygon.pAttached = 1;
    annotation.object(i).polygon.polyType = 'ground';
  elseif s>0
    % Histogram object classes in this image:
    h = hist(evidence.objNdx(setdiff(find(notDeleted),i)),[1:length(objectClasses)]);
% $$$     h = hist(evidence.objNdx(setdiff([1:Nobjects],i)),[1:length(objectClasses)]);
    
    % Compute posterior multinomial:
    mult = zeros(1,Nobjects+1);
    mult(i) = log(0);
    mult(~notDeleted) = log(0); % Remove deleted objects from consideration
    mult(end) = log(Ppart.P(s,end));
% $$$     for j = setdiff(1:Nobjects,i)
    for j = setdiff(find(notDeleted),i)
      t = evidence.objNdx(j);
      if t>0
        mult(j) = -log(Ppart.alpha_a)-(Ppart.alpha_a-1)*log(evidence.relativeArea(j)+eps)+log(Ppart.alpha_r)+(Ppart.alpha_r-1)*log(evidence.relativeOverlap(i,j)+eps)+log(1-Ppart.P(s,end))+log(Ppart.P(s,t))-log(h(t));
      else
        % Polygon is not in objectClasses:
        mult(j) = log(0);
      end
    end
    mult = exp(mult-max(mult));
    mult = mult/sum(mult);

    % Assign probability of being attached:
    annotation.object(i).polygon.pAttached = sum(mult(1:Nobjects));

    % Assign probability of which object polygon is attached to:
    annotation.object(i).polygon.pAttachedObject = mult(1:Nobjects);
  
    % Determine if object is part-of:
    if any(annotation.object(i).polygon.pAttached > threshPart)
      annotation.object(i).polygon.polyType = 'part';
      [vv,pp] = max(annotation.object(i).polygon.pAttachedObject);
      annotation.object(i).partof = annotation.object(pp).id;
    else
      annotation.object(i).polygon.polyType = 'standing';
    end
  else
    % Polygon is not in objectClasses:
    annotation.object(i).polygon.pAttached = 0;
    annotation.object(i).polygon.polyType = '';
  end
end

% Get which standing object each attached object is attached to:
allID = {annotation.object(:).id};
for i = validObjects%1:length(annotation.object)
  if strcmp(annotation.object(i).polygon.polyType,'part')
    pp = find(ismember(allID,annotation.object(i).partof));
    while strcmp(annotation.object(pp).polygon.polyType,'part')
      pp = find(ismember(allID,annotation.object(pp).partof));
    end
    annotation.object(i).partof = annotation.object(pp).id;
  end
end
