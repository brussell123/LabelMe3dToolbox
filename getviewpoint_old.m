function annotation = getviewpoint_old(annotation,params)
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

display('In getviewpoint_old.m...');

%%% LOOK INTO THIS %%%
% $$$ % Testing constraints:
% $$$ P = getCameraMatrix(annotation);
% $$$ ll = cross(P(:,3),P(:,1));
% $$$ for i = 1:length(annotation.object)
% $$$   if strcmp(annotation.object(i).world3d.type,'standingplanes')
% $$$     [cx,cy] = getContactPoints(annotation.object(i).world3d);
% $$$     [cx,cy] = LH2RH(cx,cy,imageSize);
% $$$     tt = ll'*[cx'; cy'; ones(1,length(cx))];
% $$$     if any(tt>0)
% $$$       display(tt);
% $$$     end
% $$$   end
% $$$ end

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
  k = convhull(xGround,yGround);
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
  w = t*vv;
else
  f = focalLength;
  px = 0;
  py = 0;
  Cy = 170; % Initial guess of camera height
  [xx,yy] = LH2RH(ncols/2,1,[nrows ncols]);
  yy = mean([yy max(objb)]);
  K = [f 0 px; 0 f py; 0 0 1];
  t = YtoT(max([objb maxGround]),K)-0.01;
% $$$   t = YtoT(yy,K);
% $$$   t = 0;
end

tmax = min(pi/4,YtoT(max([objb maxGround]),K));
tmin = -pi/4;

t = max(min(tmax,t),tmin);

paramsIn = [Cy t];

options = optimset;
if isfield(options,'Algorithm')
  options = optimset('NonlEqnAlgorithm','lm','Jacobian','on','Display','off','Algorithm','levenberg-marquardt');
else
  options = optimset('LargeScale','off','LevenbergMarquardt','on','Jacobian','off','Display','off');
% $$$   options = optimset('LargeScale','off','LevenbergMarquardt','on','Jacobian','on','Display','off');
% $$$   options = optimset('NonlEqnAlgorithm','lm','Jacobian','on','Display','off');
% $$$   options = optimset('LevenbergMarquardt','on','Jacobian','on','Display','off');
end

sigma_N = 0.02*nrows;%5; % Noise std dev (in pixels)

F = CostFunctionLM3D(paramsIn,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N);
display(sprintf('Initial cost: %f',F*F'/length(objh)));

% $$$ keyboard;

% Optimize:
% $$$ [paramsOpt,resnorm,residual,exitflag,output] = lsqnonlin(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,[],[],options);
[paramsOpt,resnorm,residual,exitflag,output] = lsqnonlin(@(p)CostFunctionLM3D_old(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,[0 tmin],[inf tmax],options);

% $$$ [paramsOpt,fval,exitflag,output] = fsolve(@(p)CostFunctionLM3D(p,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N),paramsIn,options);

display(paramsIn);
display(paramsOpt);

F = CostFunctionLM3D(paramsOpt,objh,objb,mu_obj,sig_obj,f,px,py,sigma_N);
display(sprintf('Final cost: %f',F*F'/length(objh)));

% Set camera output:
Cy = paramsOpt(1);
t = paramsOpt(2);

N = [0 0 0; 0 0 -1; 0 1 0];
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
