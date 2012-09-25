function [F,J] = CostFunctionLM3D(p,h,b,mu,sigma,f,px,py,sigma_N)

Cy = p(1);
t = p(2);

N = length(h); % Number of objects
F = zeros(1,N);
J = zeros(N,2);
for i = 1:N
  F(i) = Cost(Cy,t,h(i),0,b(i),mu(i),sigma(i),sigma_N,f,px,py);
  J(i,:) = Jacobian(Cy,t,h(i),0,b(i),mu(i),sigma(i),sigma_N,f,px,py);
end

F(end+1) = (Cy-170)/50;
% $$$ J(end+1,:) = [

function F = Cost(Cy,t,h,b_x,b_y,mu,sigma,sigma_N,f,px,py)
F = (Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy);

function J = Jacobian(Cy,t,h,b_x,b_y,mu,sigma,sigma_N,f,px,py)
J = [-1/2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t)))/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-1/2*(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2/(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma*(2*((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)*((sin(t)^2*px*b_y-cos(t)^2*px*b_y-py*sin(t)^2*px+py*cos(t)^2*px-2*cos(t)*f*px*sin(t))/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)+(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)^2/(f*cos(t)+py*sin(t)-sin(t)*b_y)*sin(t)-(-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)^2*(-f*sin(t)+py*cos(t)-b_y*cos(t)))+2*((-f*sin(t)+py*cos(t))/cos(t)-b_y)*((-f*cos(t)-py*sin(t))/cos(t)+(-f*sin(t)+py*cos(t))/cos(t)^2*sin(t))); h/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)-(Cy*h-(((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*mu)/((((-cos(t)*px*sin(t)*b_y+py*cos(t)*px*sin(t)+cos(t)^2*f*px+f*b_x-f*px)/cos(t)/(f*cos(t)+py*sin(t)-sin(t)*b_y)-b_x)^2+((-f*sin(t)+py*cos(t))/cos(t)-b_y)^2)^(1/2)*sigma+sigma_N*Cy)^2*sigma_N];
