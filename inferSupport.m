function sup = inferSupport(annotation,Psup,objectClasses,j)

if ~isfield(annotation,'object')
  sup = [];
  return;
end

% Parameters:
thresh_sup = 0.0001; % Threshold for support
lambda = 0.1; % Exponential distribution parameter

Nobjects = length(annotation.object);

if nargin < 4
  j = 1:Nobjects;
end

% Collect polygons:
for i = 1:Nobjects
  n = strmatch(strtrim(lower(annotation.object(i).name)), objectClasses, 'exact');
  if ~isempty(n)
    ndx(i) = n;
  else
    ndx(i) = 0;
  end
  [X{i},Y{i}] = getLMpolygon(annotation.object(i).polygon);
  [x(i),y(i)] = GetBottomPoint(X{i},Y{i});
end

Bottom = {'floor', 'ground', 'grass', 'field', 'street', 'road', ...
          'sidewalk', 'path', 'sea', 'water'};

sup = zeros(1,length(j));
P = zeros(1,length(j));
for m = 1:Nobjects
  for i = 1:length(j)
    n = j(i);
    if (m~=n) && ndx(m) && ndx(n) && ~ismember(strtrim(lower(annotation.object(n).name)),Bottom)
      pp = Psup(ndx(m),ndx(n))*exp(-lambda*SupportDist(X{m},Y{m},x(n),y(n)));
      if pp > P(i)
        P(i) = pp;
        sup(i) = m;
      end
    end
  end
end
%keyboard

% Threshold support relations:
sup(P<=thresh_sup) = 0;

% Show outputs:
[v,n] = sort(P,'descend');
for i = n
  if sup(i)
    disp(sprintf('%s supports %s (with probability = %1.2f)', ...
                 annotation.object(sup(i)).name, ...
                 annotation.object(j(i)).name,P(i)));
  end
end
