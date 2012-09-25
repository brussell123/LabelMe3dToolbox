function p = logpYCV0(v0,yc,objh,objv,mu_obj,sig_obj,w_c,mu_c,cov_c)
% p = logpYCV0(v0,yc,objh,objv,mu_obj,sig_obj,w_c,mu_c,cov_c)
%
% Log-likelihood for camera parameters.
%
% Inputs:
% v0 - Horizon line
% yc - Camera height
% objh - Object height (in pixels)
% objv - Position of foot of object (in pixels)
% mu_obj - Mean object height (in meters)
% sig_obj - Standard deviation object height (in meters)
% w_c - Prior mixture of Gaussians weights
% mu_c - Prior mixture of Gaussians means
% cov_c - Prior mixture of Gaussians covariances
%
% Outputs:
% p = -log(P(c,obj;mu, sig))

% Parameters:
noiseSigma = 0.02;

p = zeros(size(v0));    

% Derek's learned prior on camera parameters:
if nargin >= 9
  for j = 1:numel(v0)
    p(j) = 0;
    for m = 1:numel(w_c)
      p(j) = p(j) + w_c(m)*mvnpdf([v0(j) yc(j)], mu_c(m, :), cov_c(:, :, m));
    end
  end
  p = log(max(p, 1E-40));
end

Nobj = length(objh);
Npts = numel(v0);
dims = size(v0);
v0 = repmat(reshape(v0,Npts,1),1,Nobj);
yc = repmat(reshape(yc,Npts,1),1,Nobj);
objh = repmat(reshape(objh,1,Nobj),Npts,1);
objv = repmat(reshape(objv,1,Nobj),Npts,1);
mu_obj = repmat(reshape(mu_obj,1,Nobj),Npts,1);
sig_obj = repmat(reshape(sig_obj,1,Nobj),Npts,1);

mu = mu_obj.*(v0-objv)./yc;
sigma = sig_obj.*(v0-objv)./yc + noiseSigma;
p = normpdf(objh,mu,sigma);
p = sum(log(max(p,1E-40)),2);

% Prior on camera height:
p = p+log(normpdf(yc(:,1),1.7,0.5));

p = reshape(p,dims);
