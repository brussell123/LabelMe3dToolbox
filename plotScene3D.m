function plotScene3D(annotation)

clf;
range = [inf -inf inf -inf inf -inf];
for i = 1:length(annotation.object)
  if ~isempty(annotation.object(i).mesh)
    X = annotation.object(i).mesh.v3D(1,:);
    Y = annotation.object(i).mesh.v3D(2,:);
    Z = annotation.object(i).mesh.v3D(3,:);
    plot3([X X(1)],[Z Z(1)],[Y Y(1)],'b');
    hold on;
    range([1 3 5]) = min(range([1 3 5]),[min(X) min(Z) min(Y)]);
    range([2 4 6]) = max(range([2 4 6]),[max(X) max(Z) max(Y)]);
  end
end
range(3) = min(range(3),0);
axis(range);
% $$$ axis equal;
xlabel('x (cm)');
ylabel('z (cm)');
zlabel('y (cm)');
rotate3d
