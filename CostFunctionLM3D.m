function [F,J] = CostFunctionLM3D(params,xtop_in,xbot_in,mu,sigma,mu_Cy,sigma_Cy,fhat,pxhat,pyhat)

Nstanding = size(xtop_in,2);

% Camera parameters:
f = params(1);
px = params(2);
py = params(3);
wx = params(4);
wz = params(5);
Cy = params(6);
X = [reshape(params(6+[1:3*Nstanding]),3,Nstanding); ones(1,Nstanding)];

% Camera calibration:
K = [f 0 px; 0 f py; 0 0 1];

% Camera rotation:
t = sqrt(wx^2+wz^2)+eps;
nx = wx/t;
nz = wz/t;
N = [0 -nz 0; nz 0 -nx; 0 nx 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;

% Get camera center:
C = [0; Cy; 0];

% Get camera matrix:
P = K*R*[eye(3) -C];

% Compute cost vector:
x = P*X;
xtop = x(1,:)./x(3,:);
ytop = x(2,:)./x(3,:);
x = P*[X(1,:); zeros(1,Nstanding); X(3:4,:)];
xbot = x(1,:)./x(3,:);
ybot = x(2,:)./x(3,:);
F = [xtop-xtop_in(1,:); ytop-xtop_in(2,:); xbot-xbot_in(1,:); ybot-xbot_in(2,:); (X(2,:)-mu)./sigma];
F = F(:);
F(end+1) = (Cy-mu_Cy)/sigma_Cy;
F(end+1) = f-fhat;
F(end+1) = px-pxhat;
F(end+1) = py-pyhat;

if nargout>=2
  % Get derivatives of P matrix:
  [dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy] = DerivP(f,px,py,wx,wz,Cy);
  
  % Compute Jacobian matrix:
  J = ComputeJacobianStanding(X,sigma,sigma_Cy,P,dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy);
end


function J = ComputeJacobianStanding(X,sigma,sigma_Cy,P,dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy)

N = size(X,2);

Xtop = X;
Xbot = X;
Xbot(2,:) = zeros(1,N);
[Jxt,Jyt] = ComputeJxJy(Xtop,P,dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy);
[Jxb,Jyb] = ComputeJxJy(Xbot,P,dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy);
Jxb(:,8) = 0;
Jyb(:,8) = 0;

J = zeros(5*N,6+3*N);
J(:,1:6) = reshape([Jxt(:,1:6) Jyt(:,1:6) Jxb(:,1:6) Jyb(:,1:6) zeros(N,6)]',6,5*N)';
xx = reshape(1:3*N,3,N)+6;
yy = repmat(1:5:5*N,3,1);
nn = sub2ind([5*N 6+3*N],yy,xx);
J(nn) = Jxt(:,7:end)';
yy = repmat(2:5:5*N,3,1);
nn = sub2ind([5*N 6+3*N],yy,xx);
J(nn) = Jyt(:,7:end)';
yy = repmat(3:5:5*N,3,1);
nn = sub2ind([5*N 6+3*N],yy,xx);
J(nn) = Jxb(:,7:end)';
yy = repmat(4:5:5*N,3,1);
nn = sub2ind([5*N 6+3*N],yy,xx);
J(nn) = Jyb(:,7:end)';
xx = [2:3:3*N]+6;
yy = [5:5:5*N];
nn = sub2ind([5*N 6+3*N],yy,xx);
J(nn) = 1./sigma;
J = [J; zeros(1,6+3*N)];
J(end,6) = 1/sigma_Cy;
J = [J; zeros(1,6+3*N)];
J(end,1) = 1;
J = [J; zeros(1,6+3*N)];
J(end,2) = 1;
J = [J; zeros(1,6+3*N)];
J(end,3) = 1;


function [Jx,Jy] = ComputeJxJy(X,P,dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy)
N = size(X,2);

% Pre-compute:
x = P*X;
xdf = dPdf*X;
xdpx = dPdpx*X;
xdpy = dPdpy*X;
xdwx = dPdwx*X;
xdwz = dPdwz*X;
xdCy = dPdCy*X;

dxdf = (x(3,:).*xdf(1,:)-x(1,:).*xdf(3,:))./x(3,:).^2;
dydf = (x(3,:).*xdf(2,:)-x(2,:).*xdf(3,:))./x(3,:).^2;
dxdpx = (x(3,:).*xdpx(1,:)-x(1,:).*xdpx(3,:))./x(3,:).^2;
dydpx = (x(3,:).*xdpx(2,:)-x(2,:).*xdpx(3,:))./x(3,:).^2;
dxdpy = (x(3,:).*xdpy(1,:)-x(1,:).*xdpy(3,:))./x(3,:).^2;
dydpy = (x(3,:).*xdpy(2,:)-x(2,:).*xdpy(3,:))./x(3,:).^2;
dxdwx = (x(3,:).*xdwx(1,:)-x(1,:).*xdwx(3,:))./x(3,:).^2;
dydwx = (x(3,:).*xdwx(2,:)-x(2,:).*xdwx(3,:))./x(3,:).^2;
dxdwz = (x(3,:).*xdwz(1,:)-x(1,:).*xdwz(3,:))./x(3,:).^2;
dydwz = (x(3,:).*xdwz(2,:)-x(2,:).*xdwz(3,:))./x(3,:).^2;
dxdCy = (x(3,:).*xdCy(1,:)-x(1,:).*xdCy(3,:))./x(3,:).^2;
dydCy = (x(3,:).*xdCy(2,:)-x(2,:).*xdCy(3,:))./x(3,:).^2;
dxdX = (x(3,:)*P(1,1)-x(1,:)*P(3,1))./x(3,:).^2;
dydX = (x(3,:)*P(2,1)-x(2,:)*P(3,1))./x(3,:).^2;
dxdY = (x(3,:)*P(1,2)-x(1,:)*P(3,2))./x(3,:).^2;
dydY = (x(3,:)*P(2,2)-x(2,:)*P(3,2))./x(3,:).^2;
dxdZ = (x(3,:)*P(1,3)-x(1,:)*P(3,3))./x(3,:).^2;
dydZ = (x(3,:)*P(2,3)-x(2,:)*P(3,3))./x(3,:).^2;

Jx = [dxdf(:) dxdpx(:) dxdpy(:) dxdwx(:) dxdwz(:) dxdCy(:) dxdX(:) dxdY(:) dxdZ(:)];
Jy = [dydf(:) dydpx(:) dydpy(:) dydwx(:) dydwz(:) dydCy(:) dydX(:) dydY(:) dydZ(:)];


function [dPdf,dPdpx,dPdpy,dPdwx,dPdwz,dPdCy] = DerivP(f,px,py,wx,wz,Cy)

t = sqrt(wx^2+wz^2)+eps;

dPdf = [1-(1-cos(t))*wz^2/(t^2) -sin(t)*wz/t (1-cos(t))*wz/(t^2)*wx sin(t)*wz/t*Cy; ...
        sin(t)*wz/t 1-(1-cos(t))*wz^2/(t^2)-(1-cos(t))*wx^2/(t^2) -sin(t)*wx/t -(1-(1-cos(t))*wz^2/(t^2)-(1-cos(t))*wx^2/(t^2))*Cy; ...
        0 0 0 0];

dPdpx = [(1-cos(t))*wz/(t^2)*wx sin(t)*wx/t 1-(1-cos(t))*wx^2/(t^2) -sin(t)*wx/t*Cy; ...
         0 0 0 0; ...
         0 0 0 0];

dPdpy = [0 0 0 0; ...
         (1-cos(t))*wz/(t^2)*wx sin(t)*wx/t 1-(1-cos(t))*wx^2/(t^2) -sin(t)*wx/t*Cy; ...
         0 0 0 0];

dPdwx = [f*(-sin(t)*wx/(t^2)^(3/2)*wz^2+2*(1-cos(t))*wz^2/(t^2)^2*wx)+px*sin(t)*wx^2/(t^2)^(3/2)*wz-2*px*(1-cos(t))*wz/(t^2)^2*wx^2+px*(1-cos(t))*wz/(t^2) -f*cos(t)*wx/(t^2)*wz+f*sin(t)*wz/(t^2)^(3/2)*wx+px*cos(t)*wx^2/(t^2)+px*sin(t)/t-px*sin(t)*wx^2/(t^2)^(3/2) f*sin(t)*wx^2/(t^2)^(3/2)*wz-2*f*(1-cos(t))*wz/(t^2)^2*wx^2+f*(1-cos(t))*wz/(t^2)+px*(-sin(t)*wx^3/(t^2)^(3/2)-2*(1-cos(t))*wx/(t^2)+2*(1-cos(t))*wx^3/(t^2)^2) -(-f*cos(t)*wx/(t^2)*wz+f*sin(t)*wz/(t^2)^(3/2)*wx+px*cos(t)*wx^2/(t^2)+px*sin(t)/t-px*sin(t)*wx^2/(t^2)^(3/2))*Cy; ...
         f*cos(t)*wx/(t^2)*wz-f*sin(t)*wz/(t^2)^(3/2)*wx+py*sin(t)*wx^2/(t^2)^(3/2)*wz-2*py*(1-cos(t))*wz/(t^2)^2*wx^2+py*(1-cos(t))*wz/(t^2) f*(-sin(t)*wx/(t^2)^(3/2)*wz^2+2*(1-cos(t))*wz^2/(t^2)^2*wx-sin(t)*wx^3/(t^2)^(3/2)-2*(1-cos(t))*wx/(t^2)+2*(1-cos(t))*wx^3/(t^2)^2)+py*cos(t)*wx^2/(t^2)+py*sin(t)/t-py*sin(t)*wx^2/(t^2)^(3/2) -f*cos(t)*wx^2/(t^2)-f*sin(t)/t+f*sin(t)*wx^2/(t^2)^(3/2)+py*(-sin(t)*wx^3/(t^2)^(3/2)-2*(1-cos(t))*wx/(t^2)+2*(1-cos(t))*wx^3/(t^2)^2) -(f*(-sin(t)*wx/(t^2)^(3/2)*wz^2+2*(1-cos(t))*wz^2/(t^2)^2*wx-sin(t)*wx^3/(t^2)^(3/2)-2*(1-cos(t))*wx/(t^2)+2*(1-cos(t))*wx^3/(t^2)^2)+py*cos(t)*wx^2/(t^2)+py*sin(t)/t-py*sin(t)*wx^2/(t^2)^(3/2))*Cy; ...
         sin(t)*wx^2/(t^2)^(3/2)*wz-2*(1-cos(t))*wz/(t^2)^2*wx^2+(1-cos(t))*wz/(t^2) cos(t)*wx^2/(t^2)+sin(t)/t-sin(t)*wx^2/(t^2)^(3/2) -sin(t)*wx^3/(t^2)^(3/2)-2*(1-cos(t))*wx/(t^2)+2*(1-cos(t))*wx^3/(t^2)^2 -cos(t)*wx^2/(t^2)*Cy-sin(t)/t*Cy+sin(t)*wx^2/(t^2)^(3/2)*Cy];


dPdwz = [f*(-sin(t)*wz^3/(t^2)^(3/2)-2*(1-cos(t))*wz/(t^2)+2*(1-cos(t))*wz^3/(t^2)^2)+px*sin(t)*wz^2/(t^2)^(3/2)*wx+px*(1-cos(t))/(t^2)*wx-2*px*(1-cos(t))*wz^2/(t^2)^2*wx -f*cos(t)*wz^2/(t^2)-f*sin(t)/t+f*sin(t)*wz^2/(t^2)^(3/2)+px*cos(t)*wz/(t^2)*wx-px*sin(t)*wx/(t^2)^(3/2)*wz f*sin(t)*wz^2/(t^2)^(3/2)*wx+f*(1-cos(t))/(t^2)*wx-2*f*(1-cos(t))*wz^2/(t^2)^2*wx+px*(-sin(t)*wx^2/(t^2)^(3/2)*wz+2*(1-cos(t))*wz/(t^2)^2*wx^2) -(-f*cos(t)*wz^2/(t^2)-f*sin(t)/t+f*sin(t)*wz^2/(t^2)^(3/2)+px*cos(t)*wz/(t^2)*wx-px*sin(t)*wx/(t^2)^(3/2)*wz)*Cy; ...
         f*cos(t)*wz^2/(t^2)+f*sin(t)/t-f*sin(t)*wz^2/(t^2)^(3/2)+py*sin(t)*wz^2/(t^2)^(3/2)*wx+py*(1-cos(t))/(t^2)*wx-2*py*(1-cos(t))*wz^2/(t^2)^2*wx f*(-sin(t)*wz^3/(t^2)^(3/2)-2*(1-cos(t))*wz/(t^2)+2*(1-cos(t))*wz^3/(t^2)^2-sin(t)*wx^2/(t^2)^(3/2)*wz+2*(1-cos(t))*wz/(t^2)^2*wx^2)+py*cos(t)*wz/(t^2)*wx-py*sin(t)*wx/(t^2)^(3/2)*wz -f*cos(t)*wx/(t^2)*wz+f*sin(t)*wz/(t^2)^(3/2)*wx+py*(-sin(t)*wx^2/(t^2)^(3/2)*wz+2*(1-cos(t))*wz/(t^2)^2*wx^2) -(f*(-sin(t)*wz^3/(t^2)^(3/2)-2*(1-cos(t))*wz/(t^2)+2*(1-cos(t))*wz^3/(t^2)^2-sin(t)*wx^2/(t^2)^(3/2)*wz+2*(1-cos(t))*wz/(t^2)^2*wx^2)+py*cos(t)*wz/(t^2)*wx-py*sin(t)*wx/(t^2)^(3/2)*wz)*Cy; ...
         sin(t)*wx/(t^2)^(3/2)*wz^2+(1-cos(t))*wx/(t^2)-2*(1-cos(t))*wz^2/(t^2)^2*wx cos(t)*wz/(t^2)*wx-sin(t)*wx/(t^2)^(3/2)*wz -sin(t)*wx^2/(t^2)^(3/2)*wz+2*(1-cos(t))*wz/(t^2)^2*wx^2 -cos(t)*wz/(t^2)*wx*Cy+sin(t)*wx/(t^2)^(3/2)*Cy*wz];


dPdCy = [0 0 0 f*sin(t)*wz/t-px*sin(t)*wx/t; ...
         0 0 0 -f*(1-(1-cos(t))*wz^2/(t^2)-(1-cos(t))*wx^2/(t^2))-py*sin(t)*wx/t; ...
         0 0 0 -sin(t)*wx/t];









% $$$ function varargout = CostFunctionLM3D_old(varargin)
% $$$ % Evaluate cost function:
% $$$ %
% $$$ % F = CostFunctionLM3D(p,h,b,mu,sigma,f,px,py,sigma_N)
% $$$ %
% $$$ % Evaluate inequality constraints:
% $$$ %
% $$$ % F = CostFunctionLM3D(p)
% $$$ 
% $$$ J = [];
% $$$ switch nargin
% $$$  case 9
% $$$   F = EvaluateCostFunction(varargin{:});
% $$$   varargout{1} = F;
% $$$  case 6
% $$$   C = EvaluateConstraints(varargin{:});
% $$$   varargout{1} = C;
% $$$   varargout{2} = [];
% $$$ end
% $$$ 
% $$$ function C = EvaluateConstraints(p,f,px,py,x,y)
% $$$ 
% $$$ t = p(2);
% $$$ K = [f 0 px; 0 f py; 0 0 1];
% $$$ N = [0 0 0; 0 0 -1; 0 1 0];
% $$$ R = eye(3) + sin(t)*N + (1-cos(t))*N*N;
% $$$ KR = K*R;
% $$$ vx = KR(:,1);
% $$$ vz = KR(:,3);
% $$$ x = [x y ones(length(x),1)]';
% $$$ ll = cross(vz,vx);
% $$$ C = [t^2-pi^2/4 ll'*x];
% $$$ 
% $$$ function F = EvaluateCostFunction(p,h,b,mu,sigma,f,px,py,sigma_N)
% $$$ 
% $$$ Cy = p(1);
% $$$ t = p(2);
% $$$ 
% $$$ N = length(h); % Number of objects
% $$$ F = zeros(1,N);
% $$$ J = zeros(N,2);
% $$$ for i = 1:N
% $$$   F(i) = Cost(Cy,t,h(i),0,b(i),mu(i),sigma(i),sigma_N,f,px,py);
% $$$ % $$$   J(i,:) = Jacobian(Cy,t,h(i),0,b(i),mu(i),sigma(i),sigma_N,f,px,py);
% $$$ end
% $$$ 
% $$$ % Prior on camera height:
% $$$ F(end+1) = (Cy-170)/50;
% $$$ 
% $$$ F = F*F';
% $$$ 
% $$$ function F = Cost(Cy,t,h,b_x,b_y,mu,sigma,sigma_N,f,px,py)
% $$$ F = (Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy);
% $$$ 
% $$$ % $$$ function J = Jacobian(Cy,t,h,b_x,b_y,mu,sigma,sigma_N,f,px,py)
% $$$ % $$$ J = [-1/2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t)))/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-1/2*(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t))); h/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2*sigma_N];
