function pb = isThereAnyPointBelowEdge(X,Y,edges)

pb = isThereAnyPointBelow(X,Y);
xmid = 0.5*(X(edges(1,:))+X(edges(2,:)));
ymid = 0.5*(Y(edges(1,:))+Y(edges(2,:)));
pbmid = isThereAnyPointBelow(xmid,ymid);
pb = pb(edges(1,:))|pb(edges(2,:))|pbmid;
