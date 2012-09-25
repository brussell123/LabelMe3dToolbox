function n = detectShallowConcave(polygon,noPointBelow)

[X,Y] = getLMpolygon(polygon);
edges = polygon.edges;

if nargin < 2
% $$$   pb = isThereAnyPointBelow(X,Y);
% $$$   noPointBelow = 1-double(pb(edges(1,:))|pb(edges(2,:)));
  noPointBelow = 1-double(isThereAnyPointBelowEdge(X,Y,edges));
end

Nedges = size(edges,2);
k = 0;
n = [];
if Nedges < 3
  return;
end

for i = 2:Nedges-1
  if noPointBelow(i)
    yleft = Y(edges(1,i-1));
    y1 = Y(edges(1,i));
    y2 = Y(edges(2,i));
    yright = Y(edges(2,i+1));
    
% $$$     if (yleft>y1) && (yleft>y2) && (yright>y1) && (yright>y2)
    if (yleft>y1) && (yright>y2)
      k = k+1;
      n(k,:) = [i-1 i i+1];
    end
  end
end

if noPointBelow(1)
  yleft = Y(edges(1,end));
  y1 = Y(edges(1,1));
  y2 = Y(edges(2,1));
  yright = Y(edges(2,2));
  
% $$$   if (yleft>y1) && (yleft>y2) && (yright>y1) && (yright>y2)
  if (yleft>y1) && (yright>y2)
    k = k+1;
    n(k,:) = [Nedges 1 2];
  end
end

if noPointBelow(Nedges)
  yleft = Y(edges(1,end-1));
  y1 = Y(edges(1,end));
  y2 = Y(edges(2,end));
  yright = Y(edges(2,1));
  
% $$$   if (yleft>y1) && (yleft>y2) && (yright>y1) && (yright>y2)
  if (yleft>y1) && (yright>y2)
    k = k+1;
    n(k,:) = [Nedges-1 Nedges 1];
  end
end
