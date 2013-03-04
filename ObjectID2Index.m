function n = ObjectID2Index(annotation,ID)
% Inputs:
% annotation - LabelMe annotation structure
% ID - String or cell array of object IDs.
%
% Outputs:
% n - Integer indices into LabelMe annotation structure such that
%     annotation.object(n).id == ID

if isstr(ID)
  ID = {ID};
end

% Get list of object IDs:
objIDs = {annotation.object(:).id};

n = zeros(1,length(ID));
for i = 1:length(ID)
  n(i) = find(ismember(objIDs,ID{i}));
end
