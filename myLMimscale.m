function [annotation,img] = myLMimscale(annotation,img,scale,method)

% Store the original numerical format of the image and turn it into 'double'.
imgtype = whos('img');
img = single(img);

if nargin < 4
    method = 'nearest';
end

if scale ~= 1
  if nargout > 1
    % Image resampling:
    img = imresizefast(img, scale, method);
  end
  
  % add/modify image size field
  annotation.imagesize.nrows = num2str(size(img,1));
  annotation.imagesize.ncols = num2str(size(img,2));

% $$$ [annotation,img] = LMimscale(annotation,img,scale);
% $$$ if isnumeric(annotation.imagesize.nrows)
% $$$   annotation.imagesize.nrows = num2str(annotation.imagesize.nrows);
% $$$ end
% $$$ if isnumeric(annotation.imagesize.ncols)
% $$$   annotation.imagesize.ncols = num2str(annotation.imagesize.ncols);
% $$$ end

  if isfield(annotation, 'object')
    for i = 1:length(annotation.object)
      [x,y] = getLMpolygon(annotation.object(i).polygon);
      x = round(x*scale);
      y = round(y*scale);
      annotation.object(i).polygon.x = x(:);
      annotation.object(i).polygon.y = y(:);
    
      if isfield(annotation.object(i).polygon,'contact')
        [xc,yc] = getContactPoints(annotation.object(i).polygon);
        xc = round(xc*scale);
        yc = round(yc*scale);
        for j = 1:length(xc)
          annotation.object(i).polygon.contact(j).x = num2str(xc(j));
          annotation.object(i).polygon.contact(j).y = num2str(yc(j));
        end
      end
    end
  end 

  if isfield(annotation,'camera') 
    if isfield(annotation.camera,'Hy')
      Hy = str2num(annotation.camera.Hy);
      annotation.camera.Hy = num2str(Hy*scale);
    end
    if isfield(annotation.camera,'F')
      F = str2num(annotation.camera.F);
      annotation.camera.F = num2str(F*scale);
    end
  end
end

if nargout > 1
  % return the image in its original numeric format
  img = feval(imgtype.class, img);
end

function img = imresizefast(img, scaling, method, init);

if nargin<4
    init = 0;
end

if scaling > .5
    img = imresize(img, scaling, method);
else
    c = size(img,3);
    for n = 1:c
        img(:,:,n) = conv2(img(:,:,n), [1 2 1; 2 4 2; 1 2 1]/16, 'same');
    end
    img = img(init+1:2:end, init+1:2:end, :);
    %img = convn(img, [1 2 1]/4, 'same'); 
    %img = img(:,init+1:2:end,:);
    img = imresizefast(img, 2*scaling, method, 1-init);
end
