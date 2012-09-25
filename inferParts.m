function parts = inferParts(annotation,Ppart,objectClasses)

if ~isfield(annotation,'object')
  parts = [];
  return;
end

% Parameters:
thresh_par = 0.3; % Threshold for support

Nobjects = length(annotation.object);

% Collect polygons:
for i = 1:Nobjects
  n = strmatch(strtrim(lower(annotation.object(i).name)), objectClasses, 'exact');
  if ~isempty(n)
    ndx(i) = n;
  else
    ndx(i) = 0;
  end
  [X{i},Y{i}] = getLMpolygon(annotation.object(i).polygon);
end

parts = zeros(1,Nobjects);
P = zeros(1,Nobjects);
for i = 1:Nobjects
  for j = 1:Nobjects
    if (i~=j) && ndx(i) && ndx(j)
      [A,ua,area1,area2] = PolyAreas(X{i},Y{i},X{j},Y{j});
      pp = Ppart(ndx(i),ndx(j))*A/(area2+eps);
      if pp > P(j)
        P(j) = pp;
        parts(j) = i;
      end
    end
  end
end

% Threshold support relations:
parts(P<=thresh_par) = 0;

% Show outputs:
[v,n] = sort(P,'descend');
for i = n
  if parts(i)
    disp(sprintf('%s is a part of %s (with probability = %1.2f)', ...
                 annotation.object(i).name, ...
                 annotation.object(parts(i)).name,P(i)));
  end
end
