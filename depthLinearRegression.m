function [alpha,beta] = depthLinearRegression(z,zhat)
% [alpha,beta] = depthLinearRegression(z,zhat)
% 
% Find linear transformation of system depth outputs that minimizes
% relative error criterion:
%
% min  sum(abs(a*z+b-zhat)./zhat)
% a,b
%
% Inputs:
% z - System depths
% zhat - Corresponding ground truth depths
%
% Outputs:
% a,b - Linear parameters

% For now, we minimize squared relative error.  Need to minimize absolute
% relative error.
N = length(z);
alpha = (sum(z)*sum(1./zhat)+N*sum(z))/(sum(z.^2./zhat)*sum(1./zhat)+sum(z)*sum(z./zhat));
beta = (alpha*sum(z./zhat)-N)/sum(1./zhat);
