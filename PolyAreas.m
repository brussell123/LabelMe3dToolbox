function [A,ua,area1,area2] = PolyAreas(X1,Y1,X2,Y2)
% Compute intersection and union areas, as well as areas of individual
% polygons
%
% [Aint,Aunion,A1,A2] = PolyAreas(X1,Y1,X2,Y2)

  max_res = 100;
  
  min_x = min([X1(:); X2(:)]);
  max_x = max([X1(:); X2(:)]);
  min_y = min([Y1(:); Y2(:)]);
  max_y = max([Y1(:); Y2(:)]);
  
  X1 = X1-min_x;
  X2 = X2-min_x;
  Y1 = Y1-min_y;
  Y2 = Y2-min_y;
  max_x = max_x - min_x;
  max_y = max_y - min_y;
  
  max_dim = max(max_x,max_y);
  
  X1 = X1*max_res/max_dim+1; 
  X2 = X2*max_res/max_dim+1;
  Y1 = Y1*max_res/max_dim+1; 
  Y2 = Y2*max_res/max_dim+1;
  
  max_x = ceil(max_x+1); 
  max_y = ceil(max_y+1);
  
  M1 = poly2mask(X1,Y1,max_res,max_res);
  M2 = poly2mask(X2,Y2,max_res,max_res);
  
  A = sum(sum(double(M1&M2)));
  ua = sum(sum(double((M1+M2)>0)));
  area1 = sum(double(M1(:)));
  area2 = sum(double(M2(:)));
