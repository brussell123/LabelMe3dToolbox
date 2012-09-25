function [k,bb,mu_obj,sig_obj] = getObjectsWithHeightDistributions(annotation,derek,img)
% [k,bb,mu_obj,sig_obj] = getObjectsWithHeightDistributions(annotation,derek,img)
%
% Extract objects with valid height distributions.
%
% Inputs:
% annotation - LabelMe annotation structure.
% derek - Structure with object height distributions
% img - Input image
%
% Outputs:
% k - Object indices with valid height distributions
% bb - Bounding boxes of computed height extent [xmin xmax ytop ybot ycontact]
% mu_obj - Mean real-world object height.
% sig_obj - Standard deviation real-world object height.

% Extract needed information from Derek Hoiem's pre-computed structure:
objNames = derek.objNames; % Object names

if nargin < 3
  nrows = annotation.imagesize.nrows;
  ncols = annotation.imagesize.ncols;
  if isstr(nrows)
    nrows = str2num(nrows);
    ncols = str2num(ncols);
  end
else
  nrows = size(img,1);
  ncols = size(img,2);
end

Nobjects = length(annotation.object);

notDeleted = find(~isdeleted(annotation))';

% Construct object structure:
k = zeros(1,Nobjects);
bb = zeros(5,Nobjects);
mu_obj = zeros(1,Nobjects);
sig_obj = zeros(1,Nobjects);
for j = notDeleted%1:Nobjects
  % Determine whether object has computed height distribution:
  n = strmatch(strtrim(lower(annotation.object(j).name)),strtrim(lower(objNames)),'exact');

  if ~isempty(n) && (derek.mu_obj(n)~=0)
% $$$     switch annotation.object(j).polygon.polyType
    switch annotation.object(j).world3d.type
     case 'standingplanes'%'standing'
      % Get contact points:
      [cx,cy] = getContactPoints(annotation.object(j).world3d);
% $$$       [cx,cy] = getContactPoints(annotation.object(j).polygon);
      if ~isempty(cx)
        [x,y] = getLMpolygon(annotation.object(j).polygon);
        if isCleanObject(x,y,[nrows ncols])
          display(annotation.object(j).name);
          x1 = min(x); x2 = max(x); ytop = min(y);  %y2 = max(y);           
          ybot = max(cy);
          ycontact = ybot;
          k(j) = 1;
          bb(:,j) = [x1 x2 ytop ybot ycontact]';
          mu_obj(j) = derek.mu_obj(n);
          sig_obj(j) = derek.sig_obj(n);
        end
      end

% $$$      case 'part'
% $$$       % Get object it is attached to:
% $$$       nParent = getPartOfParentStanding(annotation,j);
% $$$       nParent = nParent(end);
% $$$       
% $$$       if strcmp(annotation.object(nParent).polygon.polyType,'standing')
% $$$         % Get contact points of parent object:
% $$$         cpts = find(ismember({annotation.object(nParent).polygon.pt(:).edgeType},'c'));
% $$$         
% $$$         if ~isempty(cpts)
% $$$           [x,y] = getLMpolygon(annotation.object(j).polygon);
% $$$           if isCleanObject(x,y,[nrows ncols])
% $$$             display(sprintf('%s->%s',annotation.object(j).name,annotation.object(nParent).name));
% $$$ 
% $$$             % Get extent of part:
% $$$             x1 = min(x); x2 = max(x); ytop = min(y); %[ybot,jbot] = max(y);
% $$$ % $$$             xbot = x(jbot);
% $$$             
% $$$             % Get contact point corresponding to lowest point on part:
% $$$             [xPar,yPar] = getLMpolygon(annotation.object(nParent).polygon);
% $$$ 
% $$$             %% Get bottom contact point and use this instead:
% $$$             xAll = [x1:x2];
% $$$             ycontactAll = getFootprint(xPar(cpts),yPar(cpts),xAll);
% $$$             [ycontact,ndxContact] = max(ycontactAll);
% $$$             xbot = xAll(ndxContact);
% $$$             
% $$$ % $$$             ycontact = getFootprint(xPar(cpts),yPar(cpts),xbot);
% $$$             
% $$$             clf;
% $$$             imshow(img);
% $$$             hold on;
% $$$             plot([xPar; xPar(1)],[yPar; yPar(1)],'r','LineWidth',4);
% $$$             plot([x; x(1)],[y; y(1)],'g','LineWidth',4);
% $$$             plot(xAll,ycontactAll,'y','LineWidth',4);
% $$$             plot(xbot,ycontact,'b.','MarkerSize',12);
% $$$             pause;
% $$$             
% $$$             % Set output:
% $$$           end
% $$$         end
% $$$       end
    
    end
  end
end
k = find(k);
bb = bb(:,k);
mu_obj = mu_obj(k);
sig_obj = sig_obj(k);
