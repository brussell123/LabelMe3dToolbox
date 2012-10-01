function annotation = getviewpoint(annotation,params)
% annotation = getviewpoint(annotation,params)
%
% Compute camera parameters.
%
% Inputs:
% annotation - LabelMe annotation structure
% params - Learned object height parameters
%
% Outputs:
% annotation - LabelMe annotation structure

DEBUG = 0;
if DEBUG
  annotation = getviewpoint_old(annotation,params);
  return;
end

nrows = str2num(annotation.imagesize.nrows);
ncols = str2num(annotation.imagesize.ncols);
Nobjects = length(annotation.object);

% Parameters:
focalLength = sqrt(nrows^2+ncols^2); % Predetermined focal length based
                                     % on length of image diagonal
                                     % (perhaps learn this later)

% Get object height information:
display('Computing camera calibration with the following objects:');
[bb,mu_obj,sig_obj] = getObjectsWithHeightDistributions(annotation,params,[nrows ncols]);

objh = (bb(3,:)-bb(4,:)); % Object pixel height
objb = bb(5,:); % Contact point

% Convert height distributions to centimeters:
mu_obj = 100*mu_obj;
sig_obj = 100*sig_obj;

display(sprintf('Number of objects: %d',length(objh)));

% Get points on the ground plane (ground plane and contact points):
notDeleted = find(~isdeleted(annotation))';
xGround = [];
yGround = [];
for i = notDeleted
  x = []; y = [];
  switch annotation.object(i).world3d.type
   case 'groundplane'
    [x,y] = getLMpolygon(annotation.object(i).polygon);
   case 'standingplanes'
    [x,y] = getContactPoints(annotation.object(i).world3d);
  end
  if ~isempty(x)
    [x,y] = LH2RH(x,y,[nrows ncols]);
    xGround = [xGround; x];
    yGround = [yGround; y];
  end
end

% Get convex hull of ground plane points:
if ~isempty(xGround)
  xGround(end+1) = max(xGround);
  xGround(end+1) = min(xGround);
  yGround(end+1) = min(yGround);
  yGround(end+1) = min(yGround);
  k = convhull(double(xGround),double(yGround));
  xGround = unique([xGround(k) yGround(k)],'rows');
  yGround = xGround(:,2);
  xGround = xGround(:,1);
end

% Get ground point closest to horizon line:
if ~isempty(xGround)
  maxGround = max(yGround);
else
  maxGround = -inf;
end

P = getCameraMatrix(annotation);

if ~isempty(P)
  [K,R,C] = decomposeP(P);
  f = mean([K(1) K(5)]);
  px = K(7);
  py = K(8);
  Cy = C(2);

  % Get rotation parameters from rotation matrix:
  vv = [1 0 -R(2,1)/(R(2,3)+eps)]';
  vv = vv/sqrt(vv'*vv);
  vvhat = [R(3,2)-R(2,3) R(1,3)-R(3,1) R(2,1)-R(1,2)]';
  t = atan2(vv'*vvhat,trace(R)-1);
% $$$   opt_constraints = CostFunctionLM3D([Cy t],f,px,py,xGround,yGround);
% $$$   
% $$$   if any(opt_constraints>0)
% $$$     t = YtoT(maxGround,K)-0.01;
% $$$   end
  if Cy<0
    Cy = 170;
  end
  w = t*vv;
  wx = w(1);
  wz = w(3);
else
  f = focalLength;
  px = 0;
  py = 0;
  Cy = 170; % Initial guess of camera height
  [xx,yy] = LH2RH(ncols/2,1,[nrows ncols]);
  yy = mean([yy max(objb)]);
  K = [f 0 px; 0 f py; 0 0 1];
  t = YtoT(maxGround,K)-0.01;
  wx = t;
  wz = 0;
end

paramsIn = double([f px py wx wz Cy]);

options = optimset;
if isfield(options,'Algorithm')
  options = optimset('Jacobian','on','Display','off','Algorithm','levenberg-marquardt');
else
% $$$   options = optimset('LargeScale','off','LevenbergMarquardt','on','Jacobian','on','Display','off');
  options = optimset('NonlEqnAlgorithm','lm','Jacobian','on','Display','off');
% $$$   options = optimset('NonlEqnAlgorithm','lm','Jacobian','on','Display','off','DerivativeCheck','on');
% $$$   options = optimset('LevenbergMarquardt','on','Jacobian','on','Display','off');
end

xbot_in = [zeros(1,length(objb)); objb];
xtop_in = xbot_in;
xtop_in(2,:) = xtop_in(2,:)+objh;
mu_Cy = 170;
sigma_Cy = 50;

% Back-project points:
K = [f 0 px; 0 f py; 0 0 1];
N = [0 0 0; 0 0 -1; 0 1 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;
C = [0; Cy; 0];
P = K*R*[eye(3) -C];
Xbot = P(:,[1 3 4])\[xbot_in; ones(1,length(objb))];
Zin = Xbot(2,:)./Xbot(3,:);
Xin = Xbot(1,:)./Xbot(3,:);
b = P(:,1)*Xin + P(:,3)*Zin + P(:,4)*ones(1,length(objb));
Yin = (b(2,:)-xtop_in(2,:).*b(3,:))./(xtop_in(2,:)*P(3,2)-P(2,2));
Xin = [Xin; Yin; Zin];
paramsIn = double([paramsIn Xin(:)']);

% Set regularization parameters:
fhat = f;
pxhat = px;
pyhat = py;

F = CostFunctionLM3D(paramsIn,xtop_in,xbot_in,mu_obj,sig_obj,mu_Cy,sigma_Cy,fhat,pxhat,pyhat);
display(sprintf('Initial cost: %f',F'*F));

[paramsOpt,fval,exitflag,output] = fsolve(@(p)CostFunctionLM3D(p,xtop_in,xbot_in,mu_obj,sig_obj,mu_Cy,sigma_Cy,fhat,pxhat,pyhat),paramsIn,options);

% $$$ options = optimset('NonlEqnAlgorithm','lm','Jacobian','on','Display','off');
% $$$ options = optimset('LargeScale','off','LevenbergMarquardt','on','Jacobian','on','Display','off');
% $$$ lb = [0 0 0 -pi/2 0 0 -inf*ones(1,3*length(objh))];
% $$$ ub = [inf 0 0 pi/2 0 inf inf*ones(1,3*length(objh))]
% $$$ tic;
% $$$ [paramsOpt,resnorm,residual,exitflag,output] = lsqnonlin(@(p)CostFunctionLM3D(p,xtop_in,xbot_in,mu_obj,sig_obj,mu_Cy,sigma_Cy,fhat,pxhat,pyhat),paramsIn,lb,ub,options);
% $$$ toc

% $$$ display(paramsIn);
% $$$ display(paramsOpt);

F = CostFunctionLM3D(paramsOpt,xtop_in,xbot_in,mu_obj,sig_obj,mu_Cy,sigma_Cy,fhat,pxhat,pyhat);
display(sprintf('Final cost: %f',F'*F));

% $$$ F = CostFunctionLM3D(paramsIn,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N);
% $$$ if isnan(F/length(objh))
% $$$   error('Cost is NaN');
% $$$ end
% $$$ display(sprintf('Initial cost: %f',F/length(objh)));
% $$$ % $$$ display(sprintf('Initial cost: %f',F*F'/length(objh)));
% $$$ 
% $$$ % Optimize:
% $$$ [paramsOpt,fval,exitflag,output] = fmincon(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,[],[],[],[],[],[],@(p)CostFunctionLM3D(p,f,px,py,xGround,yGround),options);
% $$$ % $$$ [paramsOpt,fval,exitflag,output] = fmincon(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,[],[],[],[],[0 tmin],[inf tmax],@(p)CostFunctionLM3D(p),options);
% $$$ % $$$ % $$$ [paramsOpt,resnorm,residual,exitflag,output] = lsqnonlin(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,[],[],options);
% $$$ % $$$ [paramsOpt,resnorm,residual,exitflag,output] = lsqnonlin(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,[0 tmin],[inf tmax],options);
% $$$ % $$$ 
% $$$ % $$$ % $$$ [paramsOpt,fval,exitflag,output] = fsolve(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,options);
% $$$ 
% $$$ display(paramsIn);
% $$$ display(paramsOpt);
% $$$ 
% $$$ F = CostFunctionLM3D(paramsOpt,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N);
% $$$ if isnan(F/length(objh))
% $$$   error('Cost is NaN');
% $$$ end
% $$$ display(sprintf('Final cost: %f',F/length(objh)));
% $$$ % $$$ display(sprintf('Final cost: %f',F*F'/length(objh)));

% Set camera output:
f = paramsOpt(1);
px = paramsOpt(2);
py = paramsOpt(3);
wx = paramsOpt(4);
wz = paramsOpt(5);
Cy = paramsOpt(6);

t = sqrt(wx^2+wz^2)+eps;
nx = wx/t;
nz = wz/t;
N = [0 -nz 0; nz 0 -nx; 0 nx 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;

K = [f 0 px; 0 f py; 0 0 1];

C = [0; Cy; 0];

P = K*R*[eye(3) -C];

annotation = setCameraMatrix(annotation,P,'centimeters');

return;

function t = YtoT(y,K)
Rx = K\[1; 0; 0];
Rx = Rx/sqrt(Rx'*Rx);
Rz = K\[0; y; 1];
Rz = Rz/sqrt(Rz'*Rz);
Ry = cross(Rz,Rx);
R = [Rx Ry Rz];
vv = [1 0 -R(2,1)/(R(2,3)+eps)]';
vv = vv/sqrt(vv'*vv);
vvhat = [R(3,2)-R(2,3) R(1,3)-R(3,1) R(2,1)-R(1,2)]';
t = atan2(vv'*vvhat,trace(R)-1);
w = t*vv;


% $$$ % get good initial value and approximate likelihood
% $$$ v0vals = [0.1:0.05:2];
% $$$ ycvals = [0.1:0.25:20];
% $$$ [v0vals, ycvals] = meshgrid(v0vals, ycvals);
% $$$ p = logpYCV0(v0vals,ycvals,objh,objv,mu_obj,sig_obj);
% $$$ [pval, ind] = max(p(:));
% $$$ v0 = v0vals(ind);
% $$$ yc = ycvals(ind);
% $$$ 
% $$$ % $$$ figure;
% $$$ % $$$ surf(v0vals,ycvals,p);
% $$$ % $$$ xlabel('horizon line');
% $$$ % $$$ ylabel('camera height');
% $$$ % $$$ 
% $$$ % $$$ objh
% $$$ % $$$ mu_obj.*(v0-objv)/yc
% $$$ % $$$ [v0 yc]
% $$$ 
% $$$ % solve precisely for v0 and yc
% $$$ options = optimset('Display','on');
% $$$ vals = fmincon(@(c) -logpYCV0(c(1),c(2),objh,objv,mu_obj,sig_obj),[v0 yc],[],[],[],[],[v0-0.5 yc-10.5],[v0+0.5 yc+10.5],[],options);
% $$$ 
% $$$ % $$$ objh
% $$$ % $$$ mu_obj.*(vals(1)-objv)/vals(2)
% $$$ % $$$ [vals(1) vals(2)]
% $$$ 
% $$$ v0 = (1-vals(1))*nrows;
% $$$ yc = vals(2);
% $$$ 
% $$$ % Make sure that horizon line is above all ground objects:
% $$$ Bottom = getListGroundObjects;
% $$$ 
% $$$ for i = 1:Nobjects
% $$$   if ismember(strtrim(lower(annotation.object(i).name)),Bottom)
% $$$     [X,Y] = getLMpolygon(annotation.object(i).polygon);
% $$$     v0 = min(min(Y)*.98, v0);
% $$$   end
% $$$ end
% $$$ 
% $$$ annotation.camera.Hy = num2str(v0);
% $$$ annotation.camera.CAM_H = num2str(100*yc); % convert to centimeters
% $$$ annotation.camera.F = num2str(focalLength);
