function annotation = insertGroundPlane(annotation)
% insertGroundPlane - Inserts a new object called "ground plane".

% Get indices of non-deleted objects:
notDeleted = find(~isdeleted(annotation))';

% Get last object ID:
maxID = annotation.object(end).id;

% Get horizon line:
[xh,yh] = getHorizonLine(annotation);

% Get image size:
imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];

annotation.object(end+1).name = 'ground';
annotation.object(end).deleted = '0';
annotation.object(end).verified = '0';
annotation.object(end).date = datestr(now);
annotation.object(end).id = num2str(str2num(maxID)+1);
annotation.object(end).polygon.pt = [];
annotation.object(end).polygon.pt(1).x = num2str(xh(1));
annotation.object(end).polygon.pt(1).y = num2str(yh(1)+1);
annotation.object(end).polygon.pt(2).x = num2str(xh(2));
annotation.object(end).polygon.pt(2).y = num2str(yh(2)+1);
annotation.object(end).polygon.pt(3).x = num2str(imageSize(2));
annotation.object(end).polygon.pt(3).y = num2str(imageSize(1));
annotation.object(end).polygon.pt(4).x = num2str(1);
annotation.object(end).polygon.pt(4).y = num2str(imageSize(1));

% Insert 3D information:
annotation.object(end).world3d.type = 'groundplane';
annotation.object(end).world3d.stale = '0';

[x,y] = getLMpolygon(annotation.object(end).polygon);
[x,y] = LH2RH(x,y,imageSize);
X = projectOntoGroundPlane(x,y,annotation);
annotation.object(end).world3d.polygon3d = setLMpolygon3D(X(1,:),X(2,:),X(3,:));

% Get plane parameters:
annotation.object(end).world3d.plane.pix = '0';
annotation.object(end).world3d.plane.piy = '1';
annotation.object(end).world3d.plane.piz = '0';
annotation.object(end).world3d.plane.piw = '0';

% Change all ground objects to be attached to ground plane:
for i = notDeleted
  if strcmp(annotation.object(i).world3d.type,'groundplane')
    annotation.object(i).world3d.type = 'part';
    annotation.object(i).world3d.parentid = annotation.object(end).id;
    annotation.object(i).world3d.rootid = annotation.object(end).id;
    for j = 1:length(annotation.object(i).world3d.polygon3d.pt)
      annotation.object(i).world3d.polygon3d.pt(j).added = '0';
      annotation.object(i).world3d.polygon3d.pt(j).planeindex.index = '0';
    end
    annotation.object(i).world3d = rmfield(annotation.object(i).world3d,'plane');
    annotation.object(i).world3d = orderfields(annotation.object(i).world3d,{'type','stale','parentid','rootid','polygon3d'});
  end
end

% Update rootid tags:
for i = notDeleted
  if strcmp(annotation.object(i).world3d.type,'part')
    % Get root object index:
    k = getPartOfParents(annotation,i);
    nRoot = k(end);

    % Set root ID:
    annotation.object(i).world3d.rootid = annotation.object(nRoot).id;
  end
end

return;
