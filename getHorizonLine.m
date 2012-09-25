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

  x1 = cross(cross(P(:,1),P(:,3)),cross([-(imageSize(2)-1)/2 0 1]',[-(imageSize(2)-1)/2 1 1]'));
  x2 = cross(cross(P(:,1),P(:,3)),cross([(imageSize(2)-1)/2 0 1]',[(imageSize(2)-1)/2 1 1]'));
  
  switch imageCoordsType
   case 'lh'
    [x1,y1] = RH2LH(x1(1)/x1(3),x1(2)/x1(3),imageSize);
    [x2,y2] = RH2LH(x2(1)/x2(3),x2(2)/x2(3),imageSize);
  end
  
  x = [x1 x2];
  y = [y1 y2];
end
