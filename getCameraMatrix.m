function [P,K,R,C] = getCameraMatrix(varargin)
% function [P,K,R,C] = getCameraMatrix(annotation)
% function [P,K,R,C] = getCameraMatrix(annotation,imageSize)
% function [P,K,R,C] = getCameraMatrix(Hy,CAM_H,F,imageSize)
% Gets camera matrix using recovered camera parameters.  Assumes
% principal point is in the center of the image.  Assumes left-handed
% world coordinate system, with Z pointing to the center (away from
% camera), X pointing right, and Y pointing up.  Assumes zero-skew and
% square pixels.
%
% The code computes the camera matrix from vanishing points in orthogonal
% directions.  The Z vanishing point is in the center on the horizon
% line.  The X vanishing point is at infinity, pointing to the right.
%
% Inputs:
% Hy - y-location in the image of the horizon line.
% CAM_H - Camera height (in world coordinates).
% F - Focal length.
% imageSize - Image size [nrows ncols]
%
% Outputs:
% P - Camera matrix.
% K - Camera calibration matrix.
% R - Rotation matrix.
% C - Camera center.

switch nargin
 case 1
  annotation = varargin{1};
  imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];
 case 2
  annotation = varargin{1};
  imageSize = varargin{2};
 case 4
  Hy = vargin{1};
  CAM_H = varargin{2};
  F = varargin{3};
  imageSize = varargin{4};
 otherwise
  error('Invalid input arguments');
end

if nargin <= 2
  if ~isfield(annotation,'camera')
    P = [];
    return;
  end
  
  if isfield(annotation.camera,'pmatrix') && isfield(annotation.camera.pmatrix,'p11')
    pmatrix = annotation.camera.pmatrix;
    P = zeros(3,4);
    if isstr(pmatrix.p11)
      P(1,1) = str2num(pmatrix.p11);
      P(1,2) = str2num(pmatrix.p12);
      P(1,3) = str2num(pmatrix.p13);
      P(1,4) = str2num(pmatrix.p14);
      P(2,1) = str2num(pmatrix.p21);
      P(2,2) = str2num(pmatrix.p22);
      P(2,3) = str2num(pmatrix.p23);
      P(2,4) = str2num(pmatrix.p24);
      P(3,1) = str2num(pmatrix.p31);
      P(3,2) = str2num(pmatrix.p32);
      P(3,3) = str2num(pmatrix.p33);
      P(3,4) = str2num(pmatrix.p34);
    else
      P(1,1) = pmatrix.p11;
      P(1,2) = pmatrix.p12;
      P(1,3) = pmatrix.p13;
      P(1,4) = pmatrix.p14;
      P(2,1) = pmatrix.p21;
      P(2,2) = pmatrix.p22;
      P(2,3) = pmatrix.p23;
      P(2,4) = pmatrix.p24;
      P(3,1) = pmatrix.p31;
      P(3,2) = pmatrix.p32;
      P(3,3) = pmatrix.p33;
      P(3,4) = pmatrix.p34;
    end
    [K,R,C] = decomposeP(P);
    return;
  end
  
  Hy = str2num(annotation.camera.Hy);
  CAM_H = str2num(annotation.camera.CAM_H);
  F = str2num(annotation.camera.F);
end

% Convert to RH coordinates:
[j,Hy] = LH2RH(0,Hy,imageSize);
% $$$ Hy = (imageSize(1)+1)/2-Hy;

% Assume principal point is in the middle of the image:
p = [0 0]';

% Camera calibration matrix:
K = [F 0 p(1); 0 F p(2); 0 0 1];

% X vanishing point:
vx = [1 0 0]';

% Z vanishing point (center of image on the horizon line):
vz = [p(1) Hy 1]';

% Get rotation matrix:
Rx = inv(K)*vx;
Rx = Rx/sqrt(Rx'*Rx);
Rz = inv(K)*vz;
Rz = Rz/sqrt(Rz'*Rz);
Ry = cross(Rz,Rx);
R = [Rx Ry Rz];

% Camera center:
C = [0 CAM_H 0]';

% Full camera matrix:
P = K*R*[eye(3) -C];
