function [x,y] = Project3D2D(varargin)
% Project 3D points to 2D points.
%
% Example 1:
% 
% annotation; % Annotation structure
% X = rand(3,10); % Generate 10 random 3D points
% [x,y] = Project3D2D(X,annotation);
%
% Example 2:
%
% P = getCameraMatrix(annotation);
% imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];
% [x,y] = Project3D2D(X,P,imageSize);

X = varargin{1};
switch nargin
 case 2
  annotation = varargin{2};
  P = getCameraMatrix(annotation);
  imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];
 case 3
  P = varargin{2};
  imageSize = varargin{3};
end

x = P*[X; ones(1,size(X,2))];
[x,y] = RH2LH(x(1,:)./x(3,:),x(2,:)./x(3,:),imageSize);

if nargout==1
  x = [x; y];
end
