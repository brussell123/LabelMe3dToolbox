function annotation = getviewpoint_ObjectHeights(annotation,params)
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

display('In getviewpoint_ObjectHeights.m...');

nrows = str2num(annotation.imagesize.nrows);
ncols = str2num(annotation.imagesize.ncols);
Nobjects = length(annotation.object);

% Parameters:
f = sqrt(nrows^2+ncols^2); % Predetermined focal length based
                           % on length of image diagonal
                           % (perhaps learn this later)
px = 0; py = 0; % Principal point
Cy = 170; % Initial guess of camera height (in centimeters)
sigma_N = 0.01*sqrt(nrows^2+ncols^2); % Labeling error std dev (in pixels)

% Get object height information:
display('Computing camera parameters with the following objects:');
[bb,mu_obj,sig_obj] = getObjectsWithHeightDistributions(annotation,params,[nrows ncols]);

objh = (bb(3,:)-bb(4,:)); % Object pixel height
objb = bb(5,:); % Contact point

% Convert height distributions to centimeters:
mu_obj = 100*mu_obj;
sig_obj = 100*sig_obj;

% Get any labeled horizon lines:
Hy = [];
notDeleted = find(~isdeleted(annotation))';
for j = notDeleted
  [x,y] = getLMpolygon(annotation.object(j).polygon);
  if strcmp(lower(annotation.object(j).originalname),'horizon line') && (length(x)==2)
    [x,y] = LH2RH(x,y,[nrows ncols]);
    Hy(end+1) = mean(y);
  end
end

if (length(objh)==0) && isempty(Hy)
  display(sprintf('No horizon line and number of valid objects: %d...cannot infer camera matrix',length(objh)));
  return;
end

display(sprintf('Number of objects: %d',length(objh)));

% Get points on the ground plane (ground plane and contact points):
notDeleted = find(~isdeleted(annotation))';
yGround = -inf;
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
    yGround = [yGround; y];
  end
end

% Get ground point closest to horizon line:
maxGround = max(yGround);

% Get bounds on camera pitch:
K = [f 0 px; 0 f py; 0 0 1];
tmin = -pi/4;
tmax = double(min(pi/4,YtoT(maxGround+1,K)));

% Initialize camera pitch to be within bounds:
t = mean([tmin tmax]);

% Initial parameters:
paramsIn = double([Cy t]);

% Optimization parameters:
options = optimset;
if isfield(options,'Algorithm')
  options = optimset('Jacobian','on','Algorithm','levenberg-marquardt','Display','off');
else
  options = optimset('Jacobian','on','NonlEqnAlgorithm','lm','Display','off');
end

[F,J] = CostFunctionLM3D_ObjectHeights(paramsIn,objh,objb,mu_obj,sig_obj,Hy,f,py,sigma_N);
display(sprintf('Initial cost: %f',F*F'));

% Optimize cost function:
[paramsOpt,resnorm,residual,exitflag,output] = lsqnonlin(@(p)CostFunctionLM3D_ObjectHeights(p,objh,objb,mu_obj,sig_obj,Hy,f,py,sigma_N),paramsIn,[0 tmin],[inf tmax],options);

display(paramsIn);
display(paramsOpt);

F = CostFunctionLM3D_ObjectHeights(paramsOpt,objh,objb,mu_obj,sig_obj,Hy,f,py,sigma_N);
display(sprintf('Final cost: %f',F*F'));

% Set camera output:
Cy = paramsOpt(1);
t = paramsOpt(2);

N = [0 0 0; 0 0 -1; 0 1 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;

K = [f 0 px; 0 f py; 0 0 1];

C = [0; Cy; 0];

P = K*R*[eye(3) -C];

% Convert camera matrix to be in LH coordinates:
imageSize = [nrows ncols];
P = [-1 0 (imageSize(2)+1)/2; 0 -1 (imageSize(1)+1)/2; 0 0 1]*P;

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
