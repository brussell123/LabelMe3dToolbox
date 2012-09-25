function annotation = contactPointsNew(annotation,sup)
%
% This function needs a full rewrite
%   - also consider edges: two nearby parallel edges are also a contact
%   cue.

thresh_sup = 50;
largeObjects  = getListLargeObjects;

Nobjects = length(annotation.object);

for i = 1:Nobjects
    [X{i},Y{i}] = getLMpolygon(annotation.object(i).polygon);
end

for i = 1:Nobjects
    % add default 0 contact
    for j = 1:length(annotation.object(i).polygon.pt)
         annotation.object(i).polygon.pt(j).contact = 0;
    end
    if sup(i) && ismember(strtrim(lower(annotation.object(i).name)),largeObjects)
        % Handle buildings/walls/etc:
% $$$         n = convhull(X{i},Y{i})';
        n = [1:length(X{i})];
        
        for j = n
            x = X{i}(j);
            y = Y{i}(j);
% $$$             if SupportDist(X{sup(i)},Y{sup(i)},x,y)<=thresh_sup && y>annotation.camera.Hy % IMPORTANT: contact points are bellow horizon line
            if y>annotation.camera.Hy % IMPORTANT: contact points are bellow horizon line
                annotation.object(i).polygon.pt(j).contact = 1;
            end
        end
    else
        % Case for all objects
        [v,n] = max(Y{i});
        if v>annotation.camera.Hy % IMPORTANT: contact points are bellow horizon line
            annotation.object(i).polygon.pt(n).contact = 1;
        end
    end
end

% This is a hack but covers an important concept: a contact point has to
% be the lowest in the same veritical line
for i = 1:Nobjects
    [x,y] = getLMpolygon(annotation.object(i).polygon);
    pb = isThereAnyPointBellow(x,y);
    for j = find(pb)'
        annotation.object(i).polygon.pt(j).contact = 0;
    end
end

% Another hack: find adjacent points that form lines with high absolute
% slope (i.e. the two poitns form a vertical line).  Do not consider the
% point at the top of the line.
for i = 1:Nobjects
  [x,y] = getLMpolygon(annotation.object(i).polygon);
  vl = isVerticalLine(x,y);
  for j = find(vl')
    annotation.object(i).polygon.pt(j).contact = 0;
  end
end

% $$$ % Handle occlusion boundaries.  Points along occlusion boundaries should
% $$$ % not be contact points.  We make sure to leave the lowest point.
% $$$ % TO DO: need to remove points that are close to the *boundary* of an
% $$$ % occluding object.  
% $$$ OCC_THRESH = 20;
% $$$ for i = 1:Nobjects
% $$$   [xi,yi] = getLMpolygon(annotation.object(i).polygon);
% $$$   c = find([annotation.object(i).polygon.pt(:).contact]);
% $$$   if ~isempty(c)
% $$$     [v,j] = max(yi(c));
% $$$     c = c(j);
% $$$     for j = 1:Nobjects
% $$$       if (i~=j) & (j~=sup(i))
% $$$         [xj,yj] = getLMpolygon(annotation.object(j).polygon);
% $$$         for k = 1:length(xi)
% $$$           if (k~=c) && annotation.object(i).polygon.pt(k).contact && (SupportDist(xj,yj,xi(k),yi(k))<OCC_THRESH)
% $$$             annotation.object(i).polygon.pt(k).contact = 0;
% $$$           end
% $$$         end
% $$$       end
% $$$     end
% $$$   end
% $$$ end

function showContactPoints(annotation,n)
% 
% shows the polygon and the contact points for object index n

[x,y] = getLMpolygon(annotation.object(n).polygon);
Cpts = [annotation.object(n).polygon.pt(:).contact];
j = find(Cpts);

clf;
plot([x; x(1)],[y; y(1)],'o-')
hold on
axis('ij')
plot(x(j), y(j), 'r*','LineWidth',4)
title(annotation.object(n).name);



function pb = isThereAnyPointBellow(x,y)
%
% A point can not be a contact point if there is another point in the same
% boundary that is bellow the current point and in the same vertical 
%
% this routine will return a indicator variable indicating if there is
% anything bellow it

x = [x;x(1)];
y = [y;y(1)];
np = length(x);
pb = zeros(np-1,1);

for n = 1:np-1
    yy = [];m=0;
    for i = 1:np-1
        if min(x(i),x(i+1))<=x(n)+0.001 & max(x(i),x(i+1))>x(n)-0.001 & i~=n
            % Insert point
            m = m+1;
            yy(m) = (y(i+1)-y(i))/(x(i+1)-x(i))*(x(n)-x(i))+y(i);
        end
    end
    if length(yy)>0
        pb(n) = y(n)<max(yy);
    end
end

function vl = isVerticalLine(x,y)
% This function returns an indicator vector if two adjacent points form a
% vertical line.  The higher point in the line is returned as "true".

vl = zeros(length(x),1);
x = [x;x(1)];
y = [y;y(1)];

SLOPE_THRESH = 4;
m = abs((y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1)));
vl = (m>SLOPE_THRESH)&(y(1:end-1)<y(2:end));
m = [m(end); m(1:end-1)];
vl = vl|((m>SLOPE_THRESH)&(y(1:end-1)<[y(end-1); y(1:end-2)]));

function d = dist2(x1,y1,x2,y2)
% x1,y1,x2,y2 are column vectors
d = repmat(x1.^2,1,length(x2))-2*x1*x2'+repmat(x2.^2,1,length(x1))' + ...
    repmat(y1.^2,1,length(y2))-2*y1*y2'+repmat(y2.^2,1,length(y1))';
