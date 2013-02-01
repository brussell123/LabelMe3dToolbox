function seg = plotPolyEdgeTypes(annotation,img,whichPolygons,method)
% Plot polygon and edge types over image.
%
% Polygon types:
%   standingplanes => red
%   part => yellow
%   groundplane => green
% Edge types:
%   occlusion => black
%   attached => gray
%   contact => white
%
% Inputs:
% annotation
% img
% whichPolygons - {valid | all}
% method - {vector | bmap}
%
% Outputs:
% seg - Output bitmap image with polygon and edge types overlayed.

if nargin < 4
  method = 'vector';%'bmap';
end
if nargout > 0
  method = 'bmap';
end
if nargin < 3
  whichPolygons = 'valid';
% $$$   whichPolygons = 'all';
end

% Parameters:
alpha = 0.5; % alpha map between polygon colors and image
cmap = [1 0 0; 0 1 0.7; 0.7 1 0];
attachedColor = [0.5 0.5 0.5];
occludedColor = [0 0 0];
contactColor = [1 1 1];

% Set line sizes:
switch method
 case 'bmap'
  sz = ceil(sqrt(size(img,1)^2+size(img,2)^2)/800);
  attachedLineWidth = sz;
  occludedLineWidth = 2*sz;
  contactLineWidth = 2*sz;
 case 'vector'
  attachedLineWidth = 1;
  occludedLineWidth = 2;
  contactLineWidth = 2;
end
  
% Remove deleted objects from consideration:
dd = isdeleted(annotation);
annotation.object = annotation.object(~dd);

% Get masks of different polygon types:
M = zeros(size(img,1),size(img,2),3);
[a,layersNdx] = LMsortlayers(annotation,img);
for i = layersNdx
  [X,Y] = getLMpolygon(annotation.object(i).polygon);
  m = poly2mask(double(X),double(Y),size(img,1),size(img,2));
  ndx = find(m(:));

  switch annotation.object(i).world3d.type
   case 'standingplanes'
    M(ndx) = cmap(1,1);
    M(ndx+size(img,1)*size(img,2)) = cmap(1,2);
    M(ndx+2*size(img,1)*size(img,2)) = cmap(1,3);
   case 'groundplane'
    M(ndx) = cmap(2,1);
    M(ndx+size(img,1)*size(img,2)) = cmap(2,2);
    M(ndx+2*size(img,1)*size(img,2)) = cmap(2,3);
   case 'part'
    parentNdx = getPartOfParents(annotation,i);
    parentNdx = parentNdx(end);
    if strcmp(whichPolygons,'all') || strcmp(annotation.object(parentNdx).world3d.type,'standingplanes')
      M(ndx) = cmap(3,1);
      M(ndx+size(img,1)*size(img,2)) = cmap(3,2);
      M(ndx+2*size(img,1)*size(img,2)) = cmap(3,3);
    end
   otherwise
    if strcmp(whichPolygons,'all')
      M(ndx) = cmap(1,1);
      M(ndx+size(img,1)*size(img,2)) = cmap(1,2);
      M(ndx+2*size(img,1)*size(img,2)) = cmap(1,3);
    end
  end
end

seg = uint8(alpha*255*M + (1-alpha)*repmat(mean(double(img),3),[1 1 3]));

if strcmp(method,'vector')
  imshow(seg);
  hold on;
end

% Plot edges for attached polygons:
for i = 1:length(annotation.object)
  if strcmp(annotation.object(i).world3d.type,'part')
    parentNdx = getPartOfParents(annotation,i);
    parentNdx = parentNdx(end);
    if strcmp(whichPolygons,'all') || strcmp(annotation.object(parentNdx).world3d.type,'standingplanes')
      % Get attached edges:
      edgesAttached = getEdges(annotation.object(i).polygon);

      switch method
       case 'bmap'
        seg = plotLineBmap(seg,edgesAttached,attachedColor,attachedLineWidth);
       case 'vector'
        plot(edgesAttached(:,[1 3])',edgesAttached(:,[2 4])','Color',attachedColor,'LineWidth',2*attachedLineWidth+1);
      end
    end
  end
end

% Plot edges for standing polygons:
for i = layersNdx
  if strcmp(whichPolygons,'all') || strcmp(annotation.object(i).world3d.type,'standingplanes')
    edgesOccluded = getEdges(annotation.object(i).polygon);
    xc = cellfun(@str2num,{annotation.object(i).world3d.contact(:).x});
    yc = cellfun(@str2num,{annotation.object(i).world3d.contact(:).y});
    edgesContact = [xc(1:end-1); yc(1:end-1); xc(2:end); yc(2:end)]';
    
    switch method
     case 'bmap'
      if ~isempty(edgesOccluded)
        seg = plotLineBmap(seg,edgesOccluded,occludedColor,occludedLineWidth);
      end
      if ~isempty(edgesContact)
        seg = plotLineBmap(seg,edgesContact,contactColor,contactLineWidth);
      end
     case 'vector'
      if ~isempty(edgesOccluded)
        plot(edgesOccluded(:,[1 3])',edgesOccluded(:,[2 4])','Color',occludedColor,'LineWidth',2*occludedLineWidth+1);
      end
      if ~isempty(edgesContact)
        plot(edgesContact(:,[1 3])',edgesContact(:,[2 4])','Color',0.99*contactColor,'LineWidth',2*contactLineWidth+1);
      end
    end
  end
end

if strcmp(method,'bmap') && (nargout==0)
  clf;
  imshow(seg);
end
