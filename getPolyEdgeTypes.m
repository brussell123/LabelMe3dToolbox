function varargout = getPolyEdgeTypes(annotation,img,type,whichPolygons)
% Inputs:
% annotation - LabelMe3D annotation structure
% img - image
% type - Any of: 'standingplanes','groundplane','part','attached_edges','contact_edges','occluded_edges'
% whichPolygons (optional) - {'valid' | 'all'}
%
% Outputs:
% seg - Binary masks for desired type
%
% Example: get binary masks for "parts" and "occluded edges"
%
% [partMask,occludedEdgesMask] = getPolyEdgeTypes(annotation,img,{'part','occluded_edges'});

if nargin < 4
  whichPolygons = 'valid';
end

if isstr(type)
  type = {type};
end

imageSize = size(img);

% Set line sizes:
sz = ceil(sqrt(imageSize(1)^2+imageSize(2)^2)/800);
attachedLineWidth = sz;
occludedLineWidth = 2*sz;
contactLineWidth = 2*sz;
  
% Remove deleted objects from consideration:
dd = isdeleted(annotation);
annotation.object = annotation.object(~dd);

SegOrder = {'standingplanes','groundplane','part','attached_edges','contact_edges','occluded_edges'};

if ~all(ismember(type,SegOrder))
  error('One of the types is invalid.');
end

% Initialize outputs:
for i = 1:length(SegOrder)
  if ismember(SegOrder{i},type)
    seg{i} = logical(zeros(imageSize(1),imageSize(2)));
  end
end

for i = 1:length(annotation.object)
  %%% Get masks of different polygon types:
  if any(ismember({'standingplanes','groundplane','part'},type))
    [X,Y] = getLMpolygon(annotation.object(i).polygon);
    mask = poly2mask(double(X),double(Y),imageSize(1),imageSize(2));
    
    switch annotation.object(i).world3d.type
     case 'standingplanes'
      if ismember('standingplanes',type)
        seg{1}(mask) = 1;
      end
     case 'groundplane'
      if ismember('groundplane',type)
        seg{2}(mask) = 1;
      end
     case 'part'
      if ismember('part',type)
        parentNdx = getPartOfParents(annotation,i);
        parentNdx = parentNdx(end);
        if strcmp(whichPolygons,'all') || strcmp(annotation.object(parentNdx).world3d.type,'standingplanes')
          seg{3}(mask) = 1;
        end
      end
     otherwise
      if strcmp(whichPolygons,'all') && ismember('standingplanes',type)
        % Render as standing plane:
        seg{1}(mask) = 1;
      end
    end
  end

  %%% Plot edges for attached polygons:
  if ismember('attached_edges',type) && strcmp(annotation.object(i).world3d.type,'part')
    parentNdx = getPartOfParents(annotation,i);
    parentNdx = parentNdx(end);
    if strcmp(whichPolygons,'all') || strcmp(annotation.object(parentNdx).world3d.type,'standingplanes')
      % Get attached edges:
      edgesAttached = getEdges(annotation.object(i).polygon);
      
      seg{4} = seg{4} | plotBinaryLine(edgesAttached,size(seg{4}),attachedLineWidth);
    end
  end

  %%% Plot edges for standing polygons:
  if strcmp(whichPolygons,'all') || strcmp(annotation.object(i).world3d.type,'standingplanes')
    edgesOccluded = getEdges(annotation.object(i).polygon);
    xc = cellfun(@str2num,{annotation.object(i).world3d.contact(:).x});
    yc = cellfun(@str2num,{annotation.object(i).world3d.contact(:).y});
    edgesContact = [xc(1:end-1); yc(1:end-1); xc(2:end); yc(2:end)]';
    
    sc = logical(zeros(imageSize(1),imageSize(2)));
    if ~isempty(edgesContact)
      sc = plotBinaryLine(edgesContact,size(sc),contactLineWidth);
    end
    if ~isempty(edgesOccluded) && ismember('occluded_edges',type)
      so = plotBinaryLine(edgesOccluded,size(seg{6}),occludedLineWidth)-sc;
      
      % Remove small drawn components:
      ll = bwlabel(so); cc = hist(ll(:),[0:max(ll(:))]);
      so = ismember(ll,find(cc(2:end)>2));

      seg{6} = seg{6} | so;
    end
    if ismember('contact_edges',type)
      seg{5} = seg{5} | sc;
    end
  end
end

% Assign outputs:
[junk,n] = ismember(type,SegOrder);
varargout = seg(n);

return;

function mask = plotBinaryLine(edges,imageSize,sz)
% Plots lines on a bitmap image.

N = ceil(max(((edges(:,3)-edges(:,1)).^2+(edges(:,4)-edges(:,2)).^2).^0.5))+10;
alpha = linspace(0,1,N);
sx = zeros(size(edges,1),N);
sy = zeros(size(edges,1),N);
for i = 1:N
  sx(:,i) = (1-alpha(i))*edges(:,1)+alpha(i)*edges(:,3);
  sy(:,i) = (1-alpha(i))*edges(:,2)+alpha(i)*edges(:,4);
end
s = round([sx(:) sy(:)]);
s = max(s,1);
s(:,1) = min(s(:,1),imageSize(2));
s(:,2) = min(s(:,2),imageSize(1));
s = unique(s,'rows');
ndx = sub2ind(imageSize,s(:,2),s(:,1));
mask = zeros(imageSize);
mask(ndx) = 1;
cB = strel('square',sz); % This is faster than disk
% $$$ cB = strel('disk',sz);
mask = imdilate(mask,cB);

return;

type = 'standingplanes';
seg = getPolyEdgeTypes(annotation,img,type);

type = 'part';
seg = getPolyEdgeTypes(annotation,img,type);

type = 'groundplane';
seg = getPolyEdgeTypes(annotation,img,type);

type = 'attached_edges';
seg = getPolyEdgeTypes(annotation,img,type);

type = 'occluded_edges';
seg = getPolyEdgeTypes(annotation,img,type);

type = 'contact_edges';
seg = getPolyEdgeTypes(annotation,img,type);

type = {'contact_edges','part'};
[ce,pa] = getPolyEdgeTypes(annotation,img,type);
