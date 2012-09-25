function annotation = getviewpoint(annotation,derek,img)
% annotation = getviewpoint(annotation,derek,img)
%
% Compute Camera Horizon/Height
%
% v0 - horizon (measured from bottom as percentage of image height)
% yc - camera height (in meters)

if isstr(annotation.imagesize.nrows)
  nrows = str2num(annotation.imagesize.nrows);
else
  nrows = annotation.imagesize.nrows;
end
if isstr(annotation.imagesize.ncols)
  ncols = str2num(annotation.imagesize.ncols);
else
  ncols = annotation.imagesize.ncols;
end

% Parameters:
focalLength = sqrt(nrows^2+ncols^2); % Predetermined focal length based
                                     % on length of image diagonal
                                     % (perhaps learn this later)
% $$$ focalLength = 800; % Predetermined focal length (perhaps learn this later)

Nobjects = length(annotation.object);

if nargin < 3
  imh = annotation.imagesize.nrows;
else
  imh = size(img,1);
  annotation.imagesize.nrows = size(img,1);
  annotation.imagesize.ncols = size(img,2);
end

% Get object height information:
display('Computing camera calibration with the following objects:');
[jObj,bb,mu_obj,sig_obj] = getObjectsWithHeightDistributions(annotation,derek);
objh = (bb(4,:)-bb(3,:))/imh; % Object pixel height
objv = 1-bb(5,:)/imh; % Contact point

display(sprintf('Number of objects: %d',length(objh)));

% get good initial value and approximate likelihood
v0vals = [0.1:0.05:2];
ycvals = [0.1:0.25:20];
[v0vals, ycvals] = meshgrid(v0vals, ycvals);
p = logpYCV0(v0vals,ycvals,objh,objv,mu_obj,sig_obj);
[pval, ind] = max(p(:));
v0 = v0vals(ind);
yc = ycvals(ind);

% $$$ figure;
% $$$ surf(v0vals,ycvals,p);
% $$$ xlabel('horizon line');
% $$$ ylabel('camera height');
% $$$ 
% $$$ objh
% $$$ mu_obj.*(v0-objv)/yc
% $$$ [v0 yc]

% solve precisely for v0 and yc
options = optimset('Display','on');
vals = fmincon(@(c) -logpYCV0(c(1),c(2),objh,objv,mu_obj,sig_obj),[v0 yc],[],[],[],[],[v0-0.5 yc-10.5],[v0+0.5 yc+10.5],[],options);

% $$$ objh
% $$$ mu_obj.*(vals(1)-objv)/vals(2)
% $$$ [vals(1) vals(2)]

v0 = (1-vals(1))*imh;
yc = vals(2);

% Make sure that horizon line is above all ground objects:
Bottom = getListGroundObjects;

for i = 1:Nobjects
  if ismember(strtrim(lower(annotation.object(i).name)),Bottom)
    [X,Y] = getLMpolygon(annotation.object(i).polygon);
    v0 = min(min(Y)*.98, v0);
  end
end

annotation.camera.Hy = num2str(v0);
annotation.camera.CAM_H = num2str(100*yc); % convert to centimeters
annotation.camera.F = num2str(focalLength);
