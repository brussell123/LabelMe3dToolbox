function polygon3D = setLMpolygon3D(x,y,z)
  
  for i = 1:length(x)
    polygon3D.pt(i).x = num2str(x(i));
    polygon3D.pt(i).y = num2str(y(i));
    polygon3D.pt(i).z = num2str(z(i));
% $$$     polygon3D.pt(i).x = num2str(round(x(i)));
% $$$     polygon3D.pt(i).y = num2str(round(y(i)));
% $$$     polygon3D.pt(i).z = num2str(round(z(i)));
  end
