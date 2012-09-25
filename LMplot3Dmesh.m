function LMplot3Dmesh(mesh)

trisurf(mesh.faces',mesh.vertices(1,:),mesh.vertices(2,:),mesh.vertices(3,:));
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
cameratoolbar('setmodeGUI','orbit');
