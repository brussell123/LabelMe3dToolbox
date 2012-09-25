function annotation = setCameraMatrix(annotation,P,units)
% Inputs:
% annotation - LabelMe annotation structure
% P - Camera matrix
% units - World units (e.g. 'centimeters', 'inches', 'none')
%
% Outputs:
% annotation

if nargin < 3
  units = 'none';
end

% Store world units:
annotation.camera.units = units;

% Store camera matrix:
annotation.camera.pmatrix.p11 = num2str(P(1,1));
annotation.camera.pmatrix.p12 = num2str(P(1,2));
annotation.camera.pmatrix.p13 = num2str(P(1,3));
annotation.camera.pmatrix.p14 = num2str(P(1,4));
annotation.camera.pmatrix.p21 = num2str(P(2,1));
annotation.camera.pmatrix.p22 = num2str(P(2,2));
annotation.camera.pmatrix.p23 = num2str(P(2,3));
annotation.camera.pmatrix.p24 = num2str(P(2,4));
annotation.camera.pmatrix.p31 = num2str(P(3,1));
annotation.camera.pmatrix.p32 = num2str(P(3,2));
annotation.camera.pmatrix.p33 = num2str(P(3,3));
annotation.camera.pmatrix.p34 = num2str(P(3,4));

% If necessary, remove old extra fields:
if isfield(annotation.camera,'Hy')
  annotation.camera = rmfield(annotation.camera,'Hy');
end
if isfield(annotation.camera,'CAM_H')
  annotation.camera = rmfield(annotation.camera,'CAM_H');
end
if isfield(annotation.camera,'F')
  annotation.camera = rmfield(annotation.camera,'F');
end
