function [X,Y,Z,valid] = Image2WorldCoords(x,y,xc,yc,P)
% [X,Y,Z] = Image2WorldCoords(x,y,xc,yc,P)
%
% This function estimates the depth of a shape that is like a folded card
% directly standing on the ground. For objects with a unique contact point,
% it assumes that the object is a flat surface parallel to the camera
% plane.
%
% Inputs:
% x,y - Image polygon coordinates
% xc,yc - Contact points
% P - 3x4 Camera matrix
%
% Outputs:
% X,Y,Z: world coordinates of the polygon
% valid - Indicates whether valid world coordinates were recovered.

N = length(x);
[K,R,C] = decomposeP(P);
valid = 1;

switch length(xc)
 case 0
  % Points living on the ground plane:
  H = P(:,[1 3 4]);
  X = H\[x(:) y(:) ones(N,1)]';
  Y = zeros(1,N);
  Z = X(2,:)./X(3,:);
  X = X(1,:)./X(3,:);
 otherwise
  if length(xc)==1
    % For objects with a unique contact point, we will assign to all the
    % points in the boundary the same depth as the contact point
    % (i.e. fronto-parallel).
    xc(2) = xc(1)+1;
    yc(2) = yc(1);
  end

% $$$   % Sort contact points from left to right:
% $$$   [xc,n] = sort(xc);
% $$$   yc = yc(n);
  
  % Get world coordinates of contact points:
  H = P(:,[1 3 4]);
  Xc = H\[xc(:) yc(:) ones(length(xc),1)]';
  Xc = [Xc(1,:)./Xc(3,:); zeros(1,length(xc)); Xc(2,:)./Xc(3,:)];
  
  % Get planes (folded cards):
  PI = zeros(4,length(xc)-1);
  Xmin = zeros(1,length(xc)-1);
  Xmax = zeros(1,length(xc)-1);
  Zmin = zeros(1,length(xc)-1);
  Zmax = zeros(1,length(xc)-1);
  for i = 1:length(xc)-1
    % Get normal vector of plane:
    n = cross(Xc(:,i+1)-Xc(:,i),[0 1 0]');
    
    % Get offset:
    pi4 = -n'*Xc(:,i);

    % Get plane parameters and bounds:
    PI(:,i) = [n; pi4];
    Xmin(i) = Xc(1,i);
    Xmax(i) = Xc(1,i+1);
    Zmin(i) = Xc(3,i);
    Zmax(i) = Xc(3,i+1);
% $$$     Xmin(i) = min(Xc(1,i:i+1));
% $$$     Xmax(i) = max(Xc(1,i:i+1));
% $$$     Zmin(i) = min(Xc(3,i:i+1));
% $$$     Zmax(i) = max(Xc(3,i:i+1));
  end
  if Xmin(1)<=Xmax(1)
    Xmin(1) = -inf;
  else
    Xmin(1) = inf;
  end
  if Xmin(end)<=Xmax(end)
    Xmax(end) = inf;
  else
    Xmax(end) = -inf;
  end
  Zmin(1) = -inf;
  Zmax(end) = inf;

  % Ray-trace image points:
  X = zeros(1,N);
  Y = zeros(1,N);
  Z = zeros(1,N);
  ismarked = logical(zeros(1,N));
  for i = 1:size(PI,2)
    lambda = -(PI(4,i)+PI(1:3,i)'*C)./(PI(1:3,i)'*R'*inv(K)*[x(:) y(:) ones(N,1)]');
    Xi = repmat(C,1,N) + repmat(lambda,3,1).*(R'*inv(K)*[x(:) y(:) ones(N,1)]');
    bb = [Xmin(i) Xmax(i)];
    n = find((Xi(1,:)-min(bb)>=-1e-6)&(Xi(1,:)-max(bb)<=1e-6));
% $$$     n = find((Xi(1,:)>=Xmin(i))&(Xi(1,:)<=Xmax(i)));
% $$$     n = find((Xi(1,:)>=Xmin(i))&(Xi(1,:)<=Xmax(i))&(Xi(3,:)>=Zmin(i))&(Xi(3,:)<=Zmax(i)));
    X(n) = Xi(1,n);
    Y(n) = Xi(2,n);
    Z(n) = Xi(3,n);
    ismarked(n) = 1;
    
% $$$     display(bb);
% $$$     display(Xi(1,foo));
% $$$     pause;
  end
  
  if ~all(ismarked)
    valid = 0;
% $$$     keyboard;
% $$$     
% $$$     XX = Xc;
% $$$     XX(2,:) = 500;%1000;
% $$$     xx = P*[XX; ones(1,size(XX,2))];
% $$$     xx = [xx(1,:)./xx(3,:); xx(2,:)./xx(3,:)];
% $$$     
% $$$     figure;
% $$$     plot(x,y);
% $$$     hold on;
% $$$     plot(x,y,'b.');
% $$$     plot(xc,yc,'rx');
% $$$     plot([xc'; xx(1,:)],[yc'; xx(2,:)],'r');
% $$$     plot(x(~ismarked),y(~ismarked),'go');
% $$$     axis ij;
% $$$     axis equal
  end
end

return;

function [X,Y,Z] = Image2WorldCoords_v01(x,y,xc,yc,imageSize,CAM_H,Hy,F)
% [X,Y,Z] = Image2WorldCoords(x,y,xc,yc,imageSize,CAM_H,Hy,F)
%
% This function estimates the depth of a shape that is like a folded card
% directly standing on the ground. For objects with a unique contact point,
% it assumes that the object is a flat surface parallel to the camera
% plane.
%
% Inputs:
% x,y - Image polygon coordinates
% xc,yc - Contact points
% imageSize - [nrows ncols] of image
% CAM_H - Camera height
% Hy - Image y coordinate of horizon line
% F - Focal length of camera
%
% Outputs:
% X,Y,Z: world coordinates of the polygon

% Parameters
MID_X = imageSize(2)/2;
MID_Y = imageSize(1)/2;

COSANG = F/(F^2+(Hy-MID_Y)^2)^0.5;
SINANG = (MID_Y-Hy)/(F^2+(Hy-MID_Y)^2)^0.5;

switch length(xc)
 case 0
  % Points living on the ground plane:
  Z = CAM_H*(F./(y-Hy)-(MID_Y-Hy)*(y-MID_Y)./(y-Hy)/F);
%  Z = CAM_H*F./(y-Hy);
 case 1
  % For objects with a unique contact point, we will assign to all the
  % points in the boundary the same depth as the contact point.
  zc = CAM_H*(F./(yc-Hy)-(MID_Y-Hy)*(yc-MID_Y)./(yc-Hy)/F);
%  zc = CAM_H*F./(yc-Hy);
  Z = zc*ones(size(x)); % constant depth
 otherwise
  % For objects with multiple contact points, use the folded card model:
  ycc = getFootprint(xc,yc,x);
  Z = CAM_H*(F./(ycc-Hy)-(MID_Y-Hy)*(ycc-MID_Y)./(ycc-Hy)/F);
%  Z = CAM_H*F./(ycc-Hy);
end

X = Z.*(x-MID_X)./(SINANG*(y-MID_Y)+COSANG*F);
%X = Z.*(x-MID_X)/F; % THIS IS WRONG!!!
Y = CAM_H-Z./(F./(y-Hy)-(MID_Y-Hy)*(y-MID_Y)./(y-Hy)/F);

if 0
  figure;
  plot3([X X(1)],[Z Z(1)],[Y Y(1)],'r.','MarkerSize',36);
  hold on;
  plot3([X X(1)],[Z Z(1)],[Y Y(1)],'b');
  xlabel('x');
  ylabel('z');
  zlabel('y');
  rotate3d
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CODE FROM groundObjectWithMesh %%%

% Parameters:
MID_X = ncols/2;
MID_Y = nrows/2;
ANG = atan2(Hy-MID_Y,F);  % angle the camera is vertically rotated
COSANG = F/(F^2+(Hy-MID_Y)^2)^0.5;
SINANG = (MID_Y-Hy)/(F^2+(Hy-MID_Y)^2)^0.5;

% Compute world coordinates for a flat horizontal surface lying on the ground 
Z = CAM_H*tan(pi/2+ANG-atan2(y-MID_Y, F));
X = Z.*(x-MID_X)./(SINANG*(y-MID_Y)+COSANG*F);
%X = Z.*(x-MID_X)/F;
Y = zeros(size(x));
TX = x;
TY = y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CODE FROM standingFoldedCardWithMesh %%%

% parameters
MID_X = ncols/2;
MID_Y = nrows/2;
ANG = atan2(Hy-MID_Y,F);  % angle the camera is vertically rotated
COSANG = F/(F^2+(Hy-MID_Y)^2)^0.5;
SINANG = (MID_Y-Hy)/(F^2+(Hy-MID_Y)^2)^0.5;

% 1) get depth for contact points with the ground, the Z coordinate are easy to compute
zc = CAM_H*tan(pi/2+ANG-atan2(yc-MID_Y, F));
% (I have copy-pasted these last two equations, it is worth checking them).


% 3) assign the other points to one of the planes defined by the edges that
% touch the ground. If a point lies outside the base, then, extrapolate the
% closest plane.
if length(j)>1
    % For objects with multiple contact points, use the folded card model:
    ycc = interp1(xc,yc,x,'linear','extrap');
    
    Z = CAM_H*tan(pi/2+ANG-atan2(ycc-MID_Y, F));
    X = Z.*(x-MID_X)./(SINANG*(y-MID_Y)+COSANG*F);
%    X = Z.*(x-MID_X)/F;
    Y = CAM_H-Z./tan(pi/2+ANG-atan2(y-MID_Y,F));
else
    % For objects with a unique contact point, we will assign to all the
    % points in the boundary the same depth as the contact point.
    Z = zc*ones(size(x)); % constant depth
    X = Z.*(x-MID_X)./(SINANG*(y-MID_Y)+COSANG*F);
%    X = Z.*(x-MID_X)/F;
    Y = CAM_H-Z./tan(pi/2+ANG-atan2(y-MID_Y,F));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ Rp = [1 0 0; 0 COSANG -SINANG; 0 SINANG COSANG];
% $$$ Sv = [1 0 0; 0 1 0; 0 0 1/F];
% $$$ f = [x-MID_X; y-MID_Y; ones(1,length(x))];
% $$$ b = inv(Sv*Rp)*f;
% $$$ a = -CAM_H./b(2,:);
% $$$ X = a.*b(1,:);
% $$$ Z = -a.*b(3,:);
% $$$ Y = zeros(size(X));
% $$$ 
% $$$ return;

