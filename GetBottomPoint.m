function [xx,yy] = GetBottomPoint(X,Y)
  [yy,low_ndx] = max(Y);
  xx = X(low_ndx);
