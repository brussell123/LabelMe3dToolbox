function thumb = LM3DpolyEdgeTypesThumbnail(annotation,img)
% Inputs:
% annotation
% img
%
% Outputs:
% thumb

% Thumbnail size (height)
tY = 96;
tY = 160;

% Scale annotation and image:
[annotation,img] = LMimscale(annotation,img,(tY+1)/size(img,1),'bilinear');

% Get polygon/edge types:
thumb3D = plotPolyEdgeTypes(annotation,img,'valid');

% sort layers
annotation = LMsortlayers(annotation, img);

if size(img,3) < 3
  img = repmat(img(:,:,1), [1 1 3]);
end

% get segmentation masks
[mask, cc, maskpol, classpol] = LMobjectmask(annotation, [size(img,1) size(img,2)]);

% Thumbnails with masks
seg = 128*ones(size(img));
if size(mask,3)>0
  M = double(colorSegments(mask, 'donotsort'))/255;
  if size(M,1)>0
    M = M + .5*repmat(sum(mask,3)==0, [1 1 3]);
    M = M / max(M(:));
    seg = M .* repmat(mean(double(img),3)/2+128, [1 1 3]);
  end
end

seg = uint8(seg);
thumb = [seg 255*ones([size(img,1),2,size(img,3)]) thumb3D];

if nargout == 0
  figure
  imshow(thumb)
end

