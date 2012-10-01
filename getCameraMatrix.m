function P = getCameraMatrix(annotation)
% function P = getCameraMatrix(annotation)
%
% Gets camera matrix for 3D scene.  
%
% Inputs:
% annotation - LabelMe3D annotation structure
%
% Outputs:
% P - 3x4 camera matrix.

P = [];

if isfield(annotation,'camera') && isfield(annotation.camera,'pmatrix') && isfield(annotation.camera.pmatrix,'p11')
  pmatrix = annotation.camera.pmatrix;
  if isstr(pmatrix.p11)
    P = [str2num(pmatrix.p11) str2num(pmatrix.p12) str2num(pmatrix.p13) str2num(pmatrix.p14); ...
         str2num(pmatrix.p21) str2num(pmatrix.p22) str2num(pmatrix.p23) str2num(pmatrix.p24); ...
         str2num(pmatrix.p31) str2num(pmatrix.p32) str2num(pmatrix.p33) str2num(pmatrix.p34)];
  else
    P = [pmatrix.p11 pmatrix.p12 pmatrix.p13 pmatrix.p14; ...
         pmatrix.p21 pmatrix.p22 pmatrix.p23 pmatrix.p24; ...
         pmatrix.p31 pmatrix.p32 pmatrix.p33 pmatrix.p34];
  end
end
