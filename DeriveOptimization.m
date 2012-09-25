% Derivation for computing cost and Jacobian for LM3D

% Camera parameters:
syms K f px py
syms C Cx Cy Cz
syms P
syms fhat lambda;
syms pxhat pyhat;
syms R t nx nz wx wz;

% Data:
syms h b_x b_y mu sigma sigma_N b

% Error vector and Jacobian:
syms F J;

% Get camera calibration matrix:
K = [f 0 px; 0 f py; 0 0 1];

% Camera rotation (assume rotation through X axis):
% $$$ t = sqrt(wx^2+wz^2);
nx = 1; %wx/t;
nz = 0; %wz/t;
N = [0 -nz 0; nz 0 -nx; 0 nx 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;

% Get camera center:
Cx = 0;
Cz = 0;
C = [Cx; Cy; Cz];

% Get camera matrix:
P = K*R*[eye(3) -C];

% Get vanishing points:
vx = P(:,1);
vy = P(:,2);
vz = P(:,3);

syms p ppx ppy
b = [b_x; b_y; 1];
p = cross(cross(vx,vz),cross(vy,b));
ppx = simplify(p(1)/p(3));
ppy = p(2)/p(3);

syms u;
u = sqrt((ppx-b_x)^2 + (ppy-b_y)^2);

F = (Cy*h-u*mu)/(u*sigma+sigma_N*Cy);

% Jacobian:
J = [diff(F,'t'); diff(F,'Cy')];
% $$$ J = [diff(F,'t'); diff(F,'Cy'); diff(F,'f'); diff(F,'px'); diff(F,'py')];

matlabFunction(F,'file','F.m');
matlabFunction(J,'file','J.m');

F = (Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy);

J = [-1/2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t)))/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-1/2*(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t))); h/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2*sigma_N];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Derivation for computing cost and Jacobian

% Camera parameters:
syms K f px py
syms C Cx Cy Cz
syms P
syms fhat lambda;
syms pxhat pyhat;
syms R t nx nz wx wz;


% Cube parameters:
syms cx cy cz dx dz ax ay az
syms sx sy sz

% Data:
syms X Y Z x y
syms p xi yi

% Error vector and Jacobian:
syms Fx Fy Jx Jy

% Get camera calibration matrix:
K = [f 0 px; 0 f py; 0 0 1];

% Camera rotation:
t = sqrt(wx^2+wz^2);
nx = wx/t;
nz = wz/t;
N = [0 -nz 0; nz 0 -nx; 0 nx 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;

% Get camera center:
C = [Cx; Cy; Cz];

% Get camera matrix:
P = K*R*[eye(3) -C];

% Get cube:
X = 0.5*[dx*sx; 2*cy*sy; dz*sz];
Ax = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
Ay = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
Az = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
A = Ax*Ay*Az;
X = A*X;
Z = X(3,:)+cz;
Y = X(2,:)+cy;
X = X(1,:)+cx;

% Projected point:
p = P*[X; Y; Z; 1];
xi = p(1)/p(3);
yi = p(2)/p(3);

% Error vector:
Fx = x-xi;
Fy = y-yi;

% Jacobian matrix:
Jx = [diff(Fx,'wx'); diff(Fx,'wz'); diff(Fx,'f'); diff(Fx,'px'); diff(Fx,'py'); diff(Fx,'cx'); diff(Fx,'cy'); diff(Fx,'cz'); diff(Fx,'dx'); diff(Fx,'dz'); diff(Fx,'ay')];
Jy = [diff(Fy,'wx'); diff(Fy,'wz'); diff(Fy,'f'); diff(Fy,'px'); diff(Fy,'py'); diff(Fy,'cx'); diff(Fy,'cy'); diff(Fy,'cz'); diff(Fy,'dx'); diff(Fy,'dz'); diff(Fy,'ay')];

matlabFunction(Fx,'file','Fx.m');
matlabFunction(Fy,'file','Fy.m');
matlabFunction(Jx,'file','Jx.m');
matlabFunction(Jy,'file','Jy.m');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below here is the derivation for computing cost and Jacobian

% Camera parameters:
syms K f px py
syms R Rx Ry Rz tx ty tz
syms C Cx Cy Cz
syms P
syms fhat lambda;
syms pxhat pyhat;


% Cube parameters:
syms cx cy cz dx dz ax ay az
% $$$ syms cx cy cz dx dy dz ax ay az
syms sx sy sz

% Data:
syms X Y Z x y
syms p xi yi

% Error vector and Jacobian:
% $$$ syms F J
syms Fx Fy Jx Jy

% Get camera calibration matrix:
K = [f 0 px; 0 f py; 0 0 1];

% Get camera rotation:
Rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
Ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
R = Rx*Ry*Rz;

% Get camera center:
C = [Cx; Cy; Cz];

% Get camera matrix:
P = K*R*[eye(3) -C];

% Get cube:
X = 0.5*[dx*sx; 2*cy*sy; dz*sz];
% $$$ X = 0.5*[dx*sx; dy*sy; dz*sz];
Ax = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
Ay = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
Az = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
A = Ax*Ay*Az;
X = A*X;
Z = X(3,:)+cz;
Y = X(2,:)+cy;
X = X(1,:)+cx;

% Projected point:
p = P*[X; Y; Z; 1];
xi = p(1)/p(3);
yi = p(2)/p(3);

% Error vector:
% $$$ F = [x-xi y-yi];
Fx = x-xi;
Fy = y-yi;

% Jacobian matrix:
Jx = [diff(Fx,'tx'); diff(Fx,'tz'); diff(Fx,'f'); diff(Fx,'px'); diff(Fx,'py'); diff(Fx,'cx'); diff(Fx,'cy'); diff(Fx,'cz'); diff(Fx,'dx'); diff(Fx,'dz'); diff(Fx,'ay')];
Jy = [diff(Fy,'tx'); diff(Fy,'tz'); diff(Fy,'f'); diff(Fy,'px'); diff(Fy,'py'); diff(Fy,'cx'); diff(Fy,'cy'); diff(Fy,'cz'); diff(Fy,'dx'); diff(Fy,'dz'); diff(Fy,'ay')];
% $$$ Jx = [diff(Fx,'tx'); diff(Fx,'tz'); diff(Fx,'f'); diff(Fx,'cx'); diff(Fx,'cy'); diff(Fx,'cz'); diff(Fx,'ay'); diff(Fx,'dx'); diff(Fx,'dz'); diff(Fx,'px'); diff(Fx,'py')];
% $$$ Jy = [diff(Fy,'tx'); diff(Fy,'tz'); diff(Fy,'f'); diff(Fy,'cx'); diff(Fy,'cy'); diff(Fy,'cz'); diff(Fy,'ay'); diff(Fy,'dx'); diff(Fy,'dz'); diff(Fy,'px'); diff(Fy,'py')];

% $$$ matlabFunction(F,'file','Cost.m');
% $$$ matlabFunction(J,'file','Jacobian.m');
matlabFunction(Fx,'file','Fx.m');
matlabFunction(Fy,'file','Fy.m');
matlabFunction(Jx,'file','Jx.m');
matlabFunction(Jy,'file','Jy.m');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below here is the derivation for computing Gradient and Hessian

% Camera parameters:
syms K f px py
syms R Rx Ry Rz tx ty tz
syms C Cx Cy Cz
syms P
syms fhat lambda;
syms pxhat pyhat;


% Cube parameters:
syms cx cy cz dx dz ax ay az
% $$$ syms cx cy cz dx dy dz ax ay az
syms sx sy sz

% Data:
syms X Y Z x y
syms p xi yi

% Cost function, gradient, and Hessian:
syms Cost Grad Hess

% Get camera calibration matrix:
K = [f 0 px; 0 f py; 0 0 1];

% Get camera rotation:
Rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
Ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
R = Rx*Ry*Rz;

% Get camera center:
C = [Cx; Cy; Cz];

% Get camera matrix:
P = K*R*[eye(3) -C];

% Get cube:
X = 0.5*[dx*sx; 2*cy*sy; dz*sz];
% $$$ X = 0.5*[dx*sx; dy*sy; dz*sz];
Ax = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
Ay = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
Az = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
A = Ax*Ay*Az;
X = A*X;
Z = X(3,:)+cz;
Y = X(2,:)+cy;
X = X(1,:)+cx;

% Projected point:
p = P*[X; Y; Z; 1];
xi = p(1)/p(3);
yi = p(2)/p(3);

% Cost function:
Cost = 0.5*(xi^2 + yi^2) - xi*x - yi*y + 0.5*lambda*(f-fhat)^2 + 0.5*lambda*(px-pxhat)^2 + 0.5*lambda*(py-pyhat)^2;

% Gradient vector:
Grad = [diff(Cost,'tx'); diff(Cost,'tz'); diff(Cost,'f'); diff(Cost,'cx'); diff(Cost,'cy'); diff(Cost,'cz'); diff(Cost,'ay'); diff(Cost,'dx'); diff(Cost,'dz'); diff(Cost,'px'); diff(Cost,'py')];
% $$$ Grad = [diff(Cost,'tx'); diff(Cost,'tz'); diff(Cost,'f'); diff(Cost,'cx'); diff(Cost,'cy'); diff(Cost,'cz'); diff(Cost,'ay'); diff(Cost,'dx'); diff(Cost,'dz')];
% $$$ Grad = [diff(Cost,'tx'); diff(Cost,'tz'); diff(Cost,'f'); diff(Cost,'cx'); diff(Cost,'cy'); diff(Cost,'cz'); diff(Cost,'ay'); diff(Cost,'dx'); diff(Cost,'dy'); diff(Cost,'dz')];
% $$$ Grad = [diff(Cost,'tx'); diff(Cost,'tz'); diff(Cost,'f')];

% Hessian matrix:
Hess = [diff(Grad,'tx') diff(Grad,'tz') diff(Grad,'f') diff(Grad,'cx') diff(Grad,'cy') diff(Grad,'cz') diff(Grad,'ay') diff(Grad,'dx') diff(Grad,'dz') diff(Grad,'px') diff(Grad,'py')];
% $$$ Hess = [diff(Grad,'tx') diff(Grad,'tz') diff(Grad,'f') diff(Grad,'cx') diff(Grad,'cy') diff(Grad,'cz') diff(Grad,'ay') diff(Grad,'dx') diff(Grad,'dz')];
% $$$ Hess = [diff(Grad,'tx') diff(Grad,'tz') diff(Grad,'f') diff(Grad,'cx') diff(Grad,'cy') diff(Grad,'cz') diff(Grad,'ay') diff(Grad,'dx') diff(Grad,'dy') diff(Grad,'dz')];
% $$$ Hess = [diff(Grad,'tx') diff(Grad,'tz') diff(Grad,'f')];


matlabFunction(Cost,'file','Cost.m');
matlabFunction(Grad,'file','Grad.m');
matlabFunction(Hess,'file','Hess.m');
