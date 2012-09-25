function theta = edgeVerticality(X,Y,edges)
% Measure verticality of edge:
% theta = 0 -> horizontal edge
% theta = 1 -> vertical edge

e1 = edges(1,:);
e2 = edges(2,:);
theta = mod(atan2(Y(e2)-Y(e1),X(e2)-X(e1)),pi);
theta = min(theta,abs(pi-theta))/0.5/pi;
