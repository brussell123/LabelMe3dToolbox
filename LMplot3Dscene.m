function LMplot3Dscene(annotation,img,varargin)
% Displays the 3D scene.
%
% The following plots 3D polygons:
%
% LMplot3Dscene(annotation)
%
% The following plots a 3D textured model:
%
% LMplot3Dscene(annotation,img)
%
% You can specify the min/max depth range:
%
% minDepth = 1; maxDepth = 10000;
% LMplot3Dscene(annotation,img,minDepth,maxDepth)

if nargin > 1
  PlotTexturedScene(annotation,img,varargin{:});
  return;
end

P = getCameraMatrix(annotation,'RH');
if isempty(P)
  display('Cannot plot since there is no camera matrix associated with this annotation file.');
  return;
end
[K,R,C] = decomposeP(P);

Ncolors = 10;
colors = hsv(Ncolors);
plot3(C(1),C(2),C(3),'ko');
hold on;
plot3(0,0,0,'ko');
plot3([C(1); 0],[C(2); 0],[C(3); 0],'k');
objStrs = [];
h = [];
for i = 1:length(annotation.object)
  if isfield(annotation.object(i).world3d,'polygon3d') && ~isempty(annotation.object(i).world3d.polygon3d)
    [X,Y,Z] = getLMpolygon3D(annotation.object(i).world3d.polygon3d);
    hi = plot3([X; X([2:end 1])]',[Y; Y([2:end 1])]',[Z; Z([2:end 1])]','Color',colors(mod(i-1,Ncolors)+1,:));
    h(end+1) = hi(1);
    objStrs{end+1} = annotation.object(i).name;
  end
end
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
legend(h,objStrs);
cameratoolbar('setmodeGUI','orbit');

function PlotTexturedScene(annotation,img,minDepth,maxDepth)

% Get depth maps:
[Xmap,Ymap,Zmap,Nmap] = getXYZmaps(annotation);
Xmap = double(Xmap);
Ymap = double(Ymap);
Zmap = double(Zmap);

% Get min/max depth to display:
MAX_DEPTH = 16777215;
zz = Zmap(Zmap~=MAX_DEPTH);
if nargin < 4
  maxDepth = prctile(zz,98);
% $$$   maxDepth = 10000;
end
if nargin < 3
  minDepth = prctile(zz,2);
% $$$   minDepth = 1;
end

% Show warped image:
% $$$ seg = plotPolyEdgeTypes(annotation3D,img,'valid');
% $$$ mask = (Zmap>minDepth).*(Zmap<maxDepth).*(mean(seg,3)>0);
mask = (Zmap>minDepth).*(Zmap<maxDepth);
Xmap2 = Xmap;
Ymap2 = Ymap;
Zmap2 = Zmap;
Xmap2(~mask)= NaN;
Ymap2(~mask)= NaN;
Zmap2(~mask)= NaN;
b = 5;
warp(Xmap2(b:end-b+1, b:end-b+1), Zmap2(b:end-b+1, b:end-b+1), Ymap2(b:end-b+1, b:end-b+1), uint8(img(b:end-b+1, b:end-b+1, : )))
axis('ij'); axis('equal')
