function LMplot3Dscene(annotation)

P = getCameraMatrix(annotation);
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
