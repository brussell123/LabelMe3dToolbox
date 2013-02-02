function [x,y] = getHorizonLine(annotation,imageCoordsType)
% Gets horizon line passing through X and Z axes for image.
%
% Inputs:
% annotation - LabelMe annotation structure.
% imageCoordsType - 'lh' or 'rh' (left-hand or right-hand coordinates)
%
% Outputs:
% x - X coordinates for horizon line
% y - Y coordinates for horizon line

if nargin < 2
  imageCoordsType = 'lh';
end

x = []; y = [];
P = getCameraMatrix(annotation);
if ~isempty(P) && isfield(annotation,'imagesize')
  imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];

% $$$   lh = cross([Project3D2D([1 0 0 0]',annotation); 1],[Project3D2D([0 0 1 0]',annotation); 1]);
  lh = cross(P(:,1),P(:,3));
  x1 = cross(lh,cross([1 1 1],[1 2 1])');
  x2 = cross(lh,cross([imageSize(2) 1 1],[imageSize(2) 2 1])');
  x = [x1(1)/x1(3) x2(1)/x2(3)];
  y = [x1(2)/x1(3) x2(2)/x2(3)];
end
