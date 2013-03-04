function annotation = AddAnglePolygonType(annotation,validObjects)
% Inputs:
% annotation - LabelMe annotation structure
% validObjects - Indices of valid objects to compute.
%
% Outputs:
% annotation - LabelMe annotation structure with added angle polygon types.

% Get indices of objects that have been labeled as an angle constraint:
nAngle = [];
for i = validObjects
  strs = regexp(lower(annotation.object(i).originalname),' ','split');
  if ismember('angle',strs)
    [x,y] = getLMpolygon(annotation.object(i).polygon);
    if length(x)==3
      nAngle(end+1) = i;
    end
  end
end

if ~isempty(nAngle)
  % Get objects that are not deleted and that are not angles:
  n = setdiff(find(isdeleted(annotation)==0)',nAngle);

  % Get ground and standing objects:
  nGround = [];
  nStanding = [];
  for i = n
    if isfield(annotation.object(i),'world3d')
      switch annotation.object(i).world3d.type
       case 'standingplanes'
        nStanding(end+1) = i;
       case 'groundplane'
        nGround(end+1) = i;
      end
    end
  end
  
  % Intersect angles with ground and standing objects:
  nStandingGround = [nStanding nGround];
  for i = nAngle
    [Xi,Yi] = getLMpolygon(annotation.object(i).polygon);
    for j = nStandingGround
      [Xj,Yj] = getLMpolygon(annotation.object(j).polygon);
      if all(inpolygon(Xi,Yi,Xj,Yj))
        % Set annotation structure with angle type:
        annotation.object(i).world3d.type = 'angle';
        annotation.object(i).world3d.stale = '0';
        annotation.object(i).world3d.angle = '90 degrees';
        annotation.object(i).world3d.rootid = annotation.object(j).id;
        if strcmp(annotation.object(j).world3d.type,'groundplane')
          annotation.object(i).world3d.planeindex = '0';
        elseif strcmp(annotation.object(j).world3d.type,'standingplanes')
          %%% NEED TO SET THIS BASED ON WHICH STANDING PLANE ANGLE
          %%% BELONGS TO:
          annotation.object(i).world3d.planeindex = '0';
        end

        break;
      end
    end
  end
end
