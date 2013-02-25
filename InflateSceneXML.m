function annotation = InflateSceneXML(annotation,notDeleted)

imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];

% Add 3D information for groundplane and standingplanes objects:
P = getCameraMatrix(annotation,'RH');
if isempty(P)
  valid = 0;
  return;
end
[K,R,C] = decomposeP(P);
for i = notDeleted%1:length(annotation.object)
  switch annotation.object(i).world3d.type
   case 'groundplane'
    % Intersect polygon points with ground plane:
    [x,y] = getLMpolygon(annotation.object(i).polygon);
    X = projectOntoGroundPlane(x,y,annotation);
    annotation.object(i).world3d.polygon3d = setLMpolygon3D(X(1,:),X(2,:),X(3,:));

    % Get plane parameters:
    annotation.object(i).world3d.plane.pix = '0';
    annotation.object(i).world3d.plane.piy = '1';
    annotation.object(i).world3d.plane.piz = '0';
    annotation.object(i).world3d.plane.piw = '0';

    % Reorder fields:
    annotation.object(i).world3d = orderfields(annotation.object(i).world3d,{'type','stale','polygon3d','plane'});
   case 'standingplanes'
    % Get planes corresponding to each contact edge by reasoning through
    % world contact lines:
    [x,y] = getLMpolygon(annotation.object(i).polygon);
    [x,y] = LH2RH(x,y,imageSize);
    xc = str2num(char({annotation.object(i).world3d.contact(:).x}));
    yc = str2num(char({annotation.object(i).world3d.contact(:).y}));
    [xc,yc] = LH2RH(xc,yc,imageSize);

    % Get plane parameters:
    H = P(:,[1 3 4]);
    if length(xc)==1
      % Make single plane frontal-parallel:
      xx = H\[xc yc 1]';
      PI = [0 0 1 -xx(2)/xx(3)]';
    else
      % Handle multiple planes:
      PI = [];
      for j = 1:length(xc)-1
        x1 = [xc(j) yc(j) 1]';
        x2 = [xc(j+1) yc(j+1) 1]';
        li = cross(x1,x2);
        lw = H'*li;
        lw = lw/sqrt(lw(1)^2+lw(2)^2);
        PI(:,j) = [lw(1) 0 lw(2) lw(3)]';
      end
    end

    % Project 2D points onto standing planes:
    [X,isAdded,nPlane,x,y] = projectOntoStandingPlanes(x,y,P,PI,xc,yc,i);
    
    if ~isempty(X)
      if any(isAdded)
        % Store 2D points:
        [x,y] = RH2LH(x,y,imageSize);
        for j = 1:length(x)
          annotation.object(i).polygon.pt(j).x = num2str(x(j));
          annotation.object(i).polygon.pt(j).y = num2str(y(j));
        end
      end
      
      % Store plane parameters:
      for j = 1:size(PI,2)
        annotation.object(i).world3d.plane(j).pix = num2str(PI(1,j));
        annotation.object(i).world3d.plane(j).piy = num2str(PI(2,j));
        annotation.object(i).world3d.plane(j).piz = num2str(PI(3,j));
        annotation.object(i).world3d.plane(j).piw = num2str(PI(4,j));
      end
    
      % Store 3D points in annotation structure:
      annotation.object(i).world3d.polygon3d = setLMpolygon3D(X(1,:),X(2,:),X(3,:));
      for j = 1:size(X,2)
        annotation.object(i).world3d.polygon3d.pt(j).added = num2str(isAdded(j));
        for k = 1:length(nPlane{j})
          annotation.object(i).world3d.polygon3d.pt(j).planeindex(k).index = num2str(nPlane{j}(k)-1);
        end
      end
      
      % Reorder fields:
      annotation.object(i).world3d = orderfields(annotation.object(i).world3d,{'type','stale','polygon3d','contact','plane'});
% $$$     else
% $$$       keyboard;
    end
   case 'part'
    % Reorder fields:
    annotation.object(i).world3d = orderfields(annotation.object(i).world3d,{'type','stale','parentid'});
   case 'none'
    % Reorder fields:
    annotation.object(i).world3d = orderfields(annotation.object(i).world3d,{'type','stale'});
  end
end

% Add 3D information for parts:
for i = notDeleted%1:length(annotation.object)
  if strcmp(annotation.object(i).world3d.type,'part')
    % Get root object index:
    k = getPartOfParents(annotation,i);
    nRoot = k(end);

    % Set root ID:
    annotation.object(i).world3d.rootid = annotation.object(nRoot).id;
    
    if isfield(annotation.object(nRoot).world3d,'polygon3d')
      % Get polygon:
      [x,y] = getLMpolygon(annotation.object(i).polygon);

      switch annotation.object(nRoot).world3d.type
       case 'standingplanes'
        % Project 2D points onto root object's standing planes:
        [x,y] = LH2RH(x,y,imageSize);
        [X,isAdded,nPlane,x,y] = projectOntoStandingPlanes(x,y,annotation,nRoot);
       case 'groundplane'
        [X,isAdded,nPlane,x,y] = projectOntoGroundPlane(x,y,annotation);
       otherwise
        error('Root object 3D type is invalid');
      end
      
      if ~isempty(X)
        if any(isAdded)
          % Store 2D points:
          [x,y] = RH2LH(x,y,imageSize);
          for j = 1:length(x)
            annotation.object(i).polygon.pt(j).x = num2str(x(j));
            annotation.object(i).polygon.pt(j).y = num2str(y(j));
          end
        end
        
        % Store 3D points in annotation structure:
        annotation.object(i).world3d.polygon3d = setLMpolygon3D(X(1,:),X(2,:),X(3,:));
        for j = 1:size(X,2)
          annotation.object(i).world3d.polygon3d.pt(j).added = num2str(isAdded(j));
          for k = 1:length(nPlane{j})
            annotation.object(i).world3d.polygon3d.pt(j).planeindex(k).index = num2str(nPlane{j}(k)-1);
          end
        end
      end
    end
  end
end
