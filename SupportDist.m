function [dist,xout,yout,alpha] = SupportDist(X,Y,xx,yy)
  
  X(end+1) = X(1);
  Y(end+1) = Y(1);
  
  %% First, test if inside polygon:
  ndx = find((yy>min(Y(2:end),Y(1:end-1))) & ...
             (yy<=max(Y(2:end),Y(1:end-1))) & ...
             (xx<=max(X(2:end),X(1:end-1))) & ...
             (Y(2:end)~=Y(1:end-1)));
  if ~isempty(ndx)
    xinter = (yy-Y(ndx)).*(X(ndx+1)-X(ndx))./(Y(ndx+1)-Y(ndx))+X(ndx);
    counter = sum(X(ndx)==X(ndx+1) | xx<= xinter);
  else
    counter = 0;
  end
  
  if mod(counter,2)
    %% We're inside of polygon; return 0 distance:
    dist = 0;
    xout = xx;
    yout = yy;
    alpha = 0;
    return;
  end
  
  %% Compute distance to inside of polygon:
  Xa = X(1:end-1); Xb = X(2:end);
  Ya = Y(1:end-1); Yb = Y(2:end);
  alpha = xx*Xa-Xa.*Xb+Xb.^2-xx*Xb+yy*Ya-Ya.*Yb+Yb.^2-yy*Yb;
  alpha = alpha./((Xa-Xb).^2+(Ya-Yb).^2+eps);
  alpha = max(min(alpha,1),0);
  
  dist = sqrt((alpha.*Xa+(1-alpha).*Xb-xx).^2 + ...
              (alpha.*Ya+(1-alpha).*Yb-yy).^2);
  [dist,min_ndx] = min(dist);
  xout = alpha(min_ndx)*Xa(min_ndx)+(1-alpha(min_ndx))*Xb(min_ndx);
  yout = alpha(min_ndx)*Ya(min_ndx)+(1-alpha(min_ndx))*Yb(min_ndx);
  alpha = alpha(min_ndx);
  
  % For support, ,let's measure just the vertical distance between
  % polygons.
  %dist = abs(max(Y)-yy);
