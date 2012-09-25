function [xx,yy,tri,j] = getMesh(x,y,npts)
%
% Input:
%   x,y: 2D polygon
%   npts: density of the grid
%
% Output:
%   xx,yy,tri: triangular mesh
%   j: indices of the boundary points in the mesh
%
% The triangular mesh can be plotted with the standard matlab command:
%   triplot(tri,xx,yy)
%
% Example:
%   x = [0 10 10 6 2]';
%   y = [0 0 5 3 6]'
%   getMesh(x,y,10);


if nargin==2
    npts = 10;
end


% Horizontal sections vertices
xm  = unique([x(:)' linspace(min(x),max(x),npts)]);

yy = []; xx = [];
for i = 1:length(xm)
    yc = intersectPolygon(x,y,xm(i));
    N = length(yc);
    if N>0
        for m = 1:N/2
            h = linspace(yc(2*m-1),yc(2*m),floor(npts/N)+2);
            yy = [yy h];
            xx = [xx xm(i)*ones(size(h))];
        end
    end
end

% make sure that all the boundary points are part of the mesh (this might
% generate triangles with area zero. So we will remove them later).
xx = [x' xx];
yy = [y' yy];

% Remove repeated points
pp = unique([xx; yy]','rows');
xx = pp(:,1)';
yy = pp(:,2)';

% Create triangular mesh
tri = delaunay(xx,yy,{'Qbb','Qc','QJ','Pp'});
% $$$ tri = delaunay(xx,yy);

% remove triangles outside the polygon
np = size(tri,1);
cx = mean(xx(tri),2);
cy = mean(yy(tri),2);
inside = inpolygon(cx,cy,x,y);
tri = tri(find(inside==1),:);

% remove triangles with zero area
areaBB = (max(xx(tri),[],2) - min(xx(tri),[],2)) .* (max(yy(tri),[],2) - min(yy(tri),[],2));
tri = tri(find(areaBB>0),:);

% recover indices boundary
j = zeros(1,length(x));
for i = 1:length(x)
    d = (xx-x(i)).^2+(yy-y(i)).^2;
    [foo, t] = min(d);
    j(i) = t(1);
end

if nargout==0
    figure
    plot(x,y,'go-')
    hold on
    triplot(tri,xx,yy)
end

return;


% $$$ function [xx,yy,tri,j] = getMesh_v02(x,y,npts)
% $$$ %
% $$$ % Input:
% $$$ %   x,y: 2D polygon
% $$$ %
% $$$ % Output:
% $$$ %   xx,yy,tri: triangular mesh
% $$$ %   j: indices of the boundary points in the mesh
% $$$ %
% $$$ % The triangular mesh can be plotted with the standard matlab command:
% $$$ %   triplot(tri,xx,yy)
% $$$ %
% $$$ % Example:
% $$$ %   x = [0 10 10 6 2]';
% $$$ %   y = [0 0 5 3 6]'
% $$$ %   getMesh(x,y);
% $$$ 
% $$$ % Horizontal sections vertices
% $$$ xm  = unique(x);
% $$$ 
% $$$ yc1 = unique(intersectPolygon(x,y,xm(1)));
% $$$ xc1 = xm(1)*ones(1,length(yc1));
% $$$ xx = xc1; yy = yc1; tri = [];
% $$$ for i = 2:length(xm)
% $$$   yc2 = unique(intersectPolygon(x,y,xm(i)));
% $$$   xc2 = xm(i)*ones(1,length(yc2));
% $$$ 
% $$$   % Triangulate planar section:
% $$$   tri = [tri; delaunay([xc1 xc2],[yc1 yc2],{'Qbb','Qc','QJ','Pp'})];
% $$$ % $$$   tri = delaunay([xc1 xc2],[yc1 yc2],{'Qt','Qbb','Qc','QJ'});
% $$$ 
% $$$   xx = [xx xc2];
% $$$   yy = [yy yc2];
% $$$   xc1 = xc2;
% $$$   yc1 = yc2;
% $$$ end
% $$$ j = 1:length(xx);
% $$$ 
% $$$ if nargout==0
% $$$     figure
% $$$     plot(x,y,'go-')
% $$$     hold on
% $$$     triplot(tri,xx,yy)
% $$$ end
% $$$ 
% $$$ 
% $$$ return;
% $$$ 
% $$$ 
