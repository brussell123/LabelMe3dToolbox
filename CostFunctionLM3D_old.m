function [F,J] = CostFunctionLM3D_old(params,h,b_y,mu,sigma,K,sigma_N)
% Inputs:
% params - Camera parameters to optimize
% h - 1xN object heights (in pixels)
% b_y - 1xN y-coordinate of contact point
% mu - 1xN mean of real-world object height (in centimeters)
% sigma - 1xN std of real-world object height (in centimeters)
% K - 3x3 camera calibration matrix
% sigma_N - std of pixel labeling error

% Get parameters to optimize:
Cy = params(1); % Camera height
t = params(2); % Camera pitch angle

% Get camera rotation matrix:
nx = 1; %wx/t;
nz = 0; %wz/t;
N = [0 -nz 0; nz 0 -nx; 0 nx 0];
R = eye(3) + sin(t)*N + (1-cos(t))*N*N;

% Get camera center:
C = [0; Cy; 0];

% Get camera matrix:
P = K*R*[eye(3) -C];

% Get horizon line:
hl = cross(P(:,1),P(:,3));

% Get foot vanishing point:
vy = P(:,2);

N = length(h); % Number of objects
F = zeros(1,N);
J = zeros(N,2);
for i = 1:N
  b_x = 0;
  b = [b_x; b_y(i); 1];
  pp = cross(hl,cross(vy,b));
  ppx = pp(1)/pp(3);
  ppy = pp(2)/pp(3);
  u = sqrt((ppx-b_x)^2 + (ppy-b_y(i))^2); % distance from horizon line to
                                          % object foot (in pixels)
  F(i) = (Cy*h(i)-u*mu(i))/(u*sigma(i)+sigma_N*Cy);
  J(i,:) = Jacobian(Cy,t,h(i),0,b_y(i),mu(i),sigma(i),sigma_N,K(1),K(7),K(8));

  % Should optimize (need to re-compute std):
  % h = (v-b)*H/Cy
end

% Regularization on camera height:
F(end+1) = (Cy-170)/50;
J(end+1,:) = [1/50 0];

return;

function J = Jacobian(Cy,t,h,b_x,b_y,mu,sigma,sigma_N,f,px,py)
J = [-1/2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t)))/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-1/2*(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t))); h/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2*sigma_N];
