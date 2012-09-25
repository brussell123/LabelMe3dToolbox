function len = edgeLength(X,Y,edges)

e1 = edges(1,:);
e2 = edges(2,:);
len = ((X(e1)-X(e2)).^2+(Y(e1)-Y(e2)).^2).^0.5;
