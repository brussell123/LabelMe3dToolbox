function annotation = removeMissedParts(annotation)
% annotation = removeMissedParts(annotation)

nrows = str2num(annotation.imagesize.nrows);
ncols = str2num(annotation.imagesize.ncols);
notDeleted = find(~isdeleted(annotation))';
Nobjects = length(annotation.object);
jobj = logical(zeros(1,Nobjects)); % indices of standing objects
for i = notDeleted%1:Nobjects
  if strcmp(annotation.object(i).polygon.polyType,'standing')
    jobj(i) = 1;
  end
end
jobj = find(jobj);
Nobj = length(jobj);

if Nobj>0
  mask = zeros(nrows,ncols,Nobj,'single');
  objarea = zeros(Nobj,1);
  for k = 1:Nobj
    [X,Y] = getLMpolygon(annotation.object(jobj(k)).polygon);
    mask(:,:,k) = poly2mask(X,Y,nrows,ncols);
  end
  intersection = reshape(mask,[nrows*ncols Nobj]);
  intersection = intersection'*intersection;
  objarea  = eps+diag(intersection);
  
  % Row i is normalized by the area of object i: intersection(i,j) =
  % proportion of object i occupied by object j
  intersection = intersection./repmat(objarea, [1 Nobj]);
  intersection = intersection .*(intersection>.05);
  
  % OnTop(i,j) = 1 => object i is on top of object j
  OnTop = intersection>0.9; % If polygon i is included in j, then obj i is on top of j

  OnTop = logical(OnTop-diag(ones(1,length(jobj))));
  [i,j] = find(OnTop);
  i = unique(jobj(i));
  if ~isempty(i)
    for k = i
      annotation.object(k).polygon.polyType = '';
    end
  end
end

