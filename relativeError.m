function err = relativeError(z,zhat)
% err = relativeError(z,zhat)
%
% Compute absolute relative error.
%
% Inputs:
% z - System depth.
% zhat - Ground truth depth.
%
% Outputs:
% err - Absolute relative error.

err = mean(abs(zhat-z)./zhat);
