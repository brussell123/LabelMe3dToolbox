function [X,isAdded,nPlane,x,y] = projectOntoStandingPlanes(varargin)

switch nargin
 case 4
  x = double(varargin{1});
  y = double(varargin{2});
  annotation = varargin{3};
  i = varargin{4};
  P = getCameraMatrix(annotation,'RH');
  imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];
  PI = getPlanes(annotation.object(i).world3d.plane);
  xc = str2num(char({annotation.object(i).world3d.contact(:).x}));
  yc = str2num(char({annotation.object(i).world3d.contact(:).y}));
  [xc,yc] = LH2RH(xc,yc,imageSize);
 case 6
  x = double(varargin{1});
  y = double(varargin{2});
  P = varargin{3};
  PI = varargin{4};
  xc = varargin{5};
  yc = varargin{6};
 case 7
  x = double(varargin{1});
  y = double(varargin{2});
  P = varargin{3};
  PI = varargin{4};
  xc = varargin{5};
  yc = varargin{6};
  i = varargin{7};
 otherwise
  error('projectOntoStandingPlanes.m: Invalid number of arguments.');
end

% Decompose camera matrix:
[K,R,C] = decomposeP(P);
Psudo = P'*inv(P*P');

% Get 3D bounds of standing planes:
[X1,X2,Z1,Z2] = standingPlaneBounds(P,PI,xc,yc);

% $$$ if i==109
% $$$   keyboard;
% $$$ end
% $$$ 
% $$$ if i==108
% $$$   keyboard;
% $$$   % point 17
% $$$ end

X = zeros(3,length(x));
isAdded = zeros(1,length(x));
nPlane = cell(1,length(x));
nPlaneOn = zeros(1,length(x));
depths = inf*ones(1,length(x));
for j = 1:size(PI,2)
  % Intersect polygon points with plane:
  xx = [x y ones(length(x),1)]';
  lambda = (-PI(:,j)'*Psudo*xx)/(PI(:,j)'*[C; 1]);
  Xlambda = Psudo*xx + repmat(lambda,4,1).*repmat([C; 1],1,length(lambda));
  Xlambda = [Xlambda(1,:)./Xlambda(4,:); Xlambda(2,:)./Xlambda(4,:); Xlambda(3,:)./Xlambda(4,:)];
  
  % Get points that lie in range of plane, is in front of camera,
  % and is the closest to the camera:
  d = sum((Xlambda-repmat(C,1,length(x))).^2,1).^0.5;
  n = find((min(X1(j),X2(j))<=Xlambda(1,:))&(max(X1(j),X2(j))>=Xlambda(1,:))&(Xlambda(3,:)>=0)&(d<depths));
% $$$   n = find((min(X1(j),X2(j))<=Xlambda(1,:))&(min(Z1(j),Z2(j))<=Xlambda(3,:))&(max(X1(j),X2(j))>=Xlambda(1,:))&(max(Z1(j),Z2(j))>=Xlambda(3,:))&(Xlambda(3,:)>=0)&(d<depths));
  
  % Assign points:
  nPlane(n) = {j};
  nPlaneOn(n) = 1;
  X(:,n) = Xlambda(:,n);
% $$$   xx = P*[Xlambda(:,n); ones(1,length(n))];
  depths(n) = d(n);
end

% $$$ figure;
% $$$ plot3(X(1,:),X(2,:),X(3,:));
% $$$ xlabel('X');
% $$$ ylabel('Y');
% $$$ zlabel('Z');
% $$$ cameratoolbar('setmodeGUI','orbit');
% $$$ 
% $$$ xx = P*[X; ones(1,size(X,2))];
% $$$ figure;
% $$$ plot(xx(1,:)./xx(3,:),xx(2,:)./xx(3,:));
% $$$ hold on;
% $$$ plot(xc,yc,'ro');

if ~any(nPlaneOn==0)
  % Refresh plane indices:
  for j = 1:length(x)
    res = PI'*[X(:,j); 1];
    n = find(abs(res)<1e-4)';
    nPlane{j} = n;
  end
  
  % Insert points between planes:
  j = 0;
  while j < size(X,2)
    j = j+1;
    k = j+1;
    if j==size(X,2)
      k = 1;
    end
    
    if all(~ismember(nPlane{j},nPlane{k}))
      % Project 3D point to image:
      xj = P*[X(:,j); 1];
      xk = P*[X(:,k); 1];
% $$$       xj = [xj(1)/xj(3); xj(2)/xj(3); 1];
% $$$       xk = [xk(1)/xk(3); xk(2)/xk(3); 1];
      
      % Get plane passing through line connecting points:
      pp = P'*cross(xj,xk);
      pp = pp/sqrt(pp(1:3)'*pp(1:3));
      
      % Get indices of adjacent planes:
      dd = abs(repmat(nPlane{j}',1,length(nPlane{k}))-repmat(nPlane{k},length(nPlane{j}),1));
      [vv,nn] = min(dd(:));
      [ji,ki] = ind2sub(size(dd),nn);
      ii = linspace(nPlane{j}(ji),nPlane{k}(ki),abs(nPlane{j}(ji)-nPlane{k}(ki))+1);

      %%%%%%%%%%%%%%%%% NEW CODE %%%%%%%%%%%%%%%%%%%%%%%%%
% $$$       if exist('i','var') && (i==13)
% $$$         keyboard;
% $$$       end
      
      % Get contact point at boundary between planes:
      nc = intersect([ii(1) ii(1)+1],[ii(2) ii(2)+1]);
      
      % Get plane passing through contact point orthogonal to first plane:
      H = P(:,[1 3 4]);
      Xc = H\[xc(nc); yc(nc); 1];
      PI_c = [-PI(3,ii(1)) 0 PI(1,ii(1)) PI(3,ii(1))*Xc(1)/Xc(3)-PI(1,ii(1))*Xc(2)/Xc(3)]';

      % Get 3D point that is intersection of three planes:
      Xp = null([pp PI(:,ii(1)) PI_c]');
      Xp = [Xp(1)/Xp(4) Xp(2)/Xp(4) Xp(3)/Xp(4)]';
      
      doDisplay = 0;
      if doDisplay
        figure;
        plot3([X(1,:) X(1,1)],[X(2,:) X(2,1)],[X(3,:) X(3,1)]);
        hold on;
        plot3(Xp(1),Xp(2),Xp(3),'r+');        
      end        
      
      % Insert point into list:
      xx = P*[Xp; 1];
      if j<k
        x = [x(1:j); xx(1)/xx(3); x(j+1:end)];
        y = [y(1:j); xx(2)/xx(3); y(j+1:end)];
        X = [X(:,1:j) Xp X(:,j+1:end)];
        nPlane = [nPlane(1:j) [ii(1) ii(2)] nPlane(j+1:end)];
        isAdded = [isAdded(1:j) 1 isAdded(j+1:end)];
      else
        x = [x(1:j); xx(1)/xx(3)];
        y = [y(1:j); xx(2)/xx(3)];
        X = [X(:,1:j) Xp];
        nPlane = [nPlane(1:j) [ii(1) ii(2)]];
        isAdded = [isAdded(1:j) 1];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
% $$$       % Check to make sure that planes are not the same:
% $$$       n1 = PI(1:3,ii(1));
% $$$       n2 = PI(1:3,ii(2));
% $$$       n1 = n1/sqrt(n1'*n1);
% $$$       n2 = n2/sqrt(n2'*n2);
% $$$       if abs(1-n1'*n2)>1e-4
% $$$         % Get 3D point that is intersection of three planes:
% $$$         Xp = null([pp PI(:,ii(1)) PI(:,ii(2))]');
% $$$         Xp = [Xp(1)/Xp(4) Xp(2)/Xp(4) Xp(3)/Xp(4)]';
% $$$         
% $$$         % Insert point into list:
% $$$         xx = P*[Xp; 1];
% $$$         if j<k
% $$$           x = [x(1:j); xx(1)/xx(3); x(j+1:end)];
% $$$           y = [y(1:j); xx(2)/xx(3); y(j+1:end)];
% $$$           X = [X(:,1:j) Xp X(:,j+1:end)];
% $$$           nPlane = [nPlane(1:j) [ii(1) ii(2)] nPlane(j+1:end)];
% $$$           isAdded = [isAdded(1:j) 1 isAdded(j+1:end)];
% $$$         else
% $$$           x = [x(1:j); xx(1)/xx(3)];
% $$$           y = [y(1:j); xx(2)/xx(3)];
% $$$           X = [X(:,1:j) Xp];
% $$$           nPlane = [nPlane(1:j) [ii(1) ii(2)]];
% $$$           isAdded = [isAdded(1:j) 1];
% $$$         end
% $$$       end
    end
  end
else
  X = [];
  isAdded = [];
  nPlane = [];
end


function [X1,X2,Z1,Z2] = standingPlaneBounds(P,PI,xc,yc)
% Gets 3D bounds of standing planes along X and Z axes.

H = P(:,[1 3 4]);
if size(PI,2)==1
  X1 = -inf; X2 = inf;
  Z1 = -inf; Z2 = inf;
else
  X1 = []; Z1 = []; X2 = []; Z2 = [];
  for j = 1:size(PI,2)
    x1 = [xc(j) yc(j) 1]';
    x2 = [xc(j+1) yc(j+1) 1]';
    xx1 = H\x1;
    X1(j) = xx1(1)/xx1(3);
    Z1(j) = xx1(2)/xx1(3);
    xx2 = H\x2;
    X2(j) = xx2(1)/xx2(3);
    Z2(j) = xx2(2)/xx2(3);
  end
  
  % Adjust first and last planes to extend to infinity
  Xmin = (X1(1)-X2(1)-eps)/0; Zmin = (Z1(1)-Z2(1)-eps)/0;
  Xmax = (X2(end)-X1(end)+eps)/0; Zmax = (Z2(end)-Z1(end)+eps)/0;
  X1(1) = Xmin; Z1(1) = Zmin;
  X2(end) = Xmax; Z2(end) = Zmax;
end
