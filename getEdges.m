function edges = getEdges(polygon,j)
% Inputs:
% polygon
% j
%
% Outputs:
% edges

[X,Y] = getLMpolygon(polygon);
allEdges = [[1:length(X)]; [2:length(X) 1]];
if nargin < 2
  j = 1:size(allEdges,2);
end
edges = [X(allEdges(1,j)) Y(allEdges(1,j)) X(allEdges(2,j)) Y(allEdges(2,j))];

return;

% $$$ allEdges = [[1:length(polygon.pt)]; [2:length(polygon.pt) 1]];
% $$$ 
% $$$ if nargin < 2
% $$$   j = 1:size(allEdges,2);
% $$$ end
% $$$ 
% $$$ [X,Y] = getLMpolygon(polygon);
% $$$ edges = [X(allEdges(1,j)) Y(allEdges(1,j)) X(allEdges(2,j)) Y(allEdges(2,j))];
