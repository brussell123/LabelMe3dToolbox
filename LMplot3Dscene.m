function LMplot3Dscene(annotation,img)

if nargin > 1
  PlotTexturedScene(annotation,img);
  return;
end

P = getCameraMatrix(annotation,'RH');
[K,R,C] = decomposeP(P);

Ncolors = 10;
colors = hsv(Ncolors);
% $$$ figure;
% $$$ plot(C(1),C(3),'ko');
plot3(C(1),C(2),C(3),'ko');
hold on;
% $$$ plot(0,0,'ko');
plot3(0,0,0,'ko');
% $$$ plot([C(1); 0],[C(3); 0],'k');
plot3([C(1); 0],[C(2); 0],[C(3); 0],'k');
objStrs = [];
h = [];
for i = 1:length(annotation.object)
  if isfield(annotation.object(i).world3d,'polygon3d') && ~isempty(annotation.object(i).world3d.polygon3d)
    [X,Y,Z] = getLMpolygon3D(annotation.object(i).world3d.polygon3d);
% $$$     hi = plot([X X([2:end 1])]',[Z Z([2:end 1])]','Color',colors(mod(i-1,Ncolors)+1,:));
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

function PlotTexturedScene(annotation,img)

[Xmap,Ymap,Zmap,Nmap] = getXYZmaps(annotation);
Xmap = double(Xmap);
Ymap = double(Ymap);
Zmap = double(Zmap);

% Show warped image:
% $$$ seg = plotPolyEdgeTypes(annotation3D,img,'valid');
% $$$ mask = (Zmap>1).*(Zmap<10000).*(mean(seg,3)>0);
mask = (Zmap>1).*(Zmap<10000);
Xmap2 = Xmap;
Ymap2 = Ymap;
Zmap2 = Zmap;
Xmap2(~mask)= NaN;
Ymap2(~mask)= NaN;
Zmap2(~mask)= NaN;
b = 5;
figure
warp(Xmap2(b:end-b+1, b:end-b+1), Zmap2(b:end-b+1, b:end-b+1), Ymap2(b:end-b+1, b:end-b+1), uint8(img(b:end-b+1, b:end-b+1, : )))
axis('ij'); axis('equal')
