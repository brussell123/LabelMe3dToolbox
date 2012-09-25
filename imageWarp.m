function imgWarp = imageWarp(Ox, Oy, imgOriginal, Wx, Wy, imgWarp, boundary)
%
% Input:
%   Ox, Oy: landmarks in input image
%   imgOriginal: original image to be warped
%   Wx, Wy: landmarks in warped image
%   imgWard: target image 
%   boundary: indices of landmarks defining the boundary of the polygon
%
% Output:
%   imgWarp


% April, 15, 2002
% A. Torralba
% torralba@ai.mit.edu

[nrows, ncols, cc] = size(imgWarp);

% coordinates in warped image
[denseWx, denseWy] = meshgrid(floor(min(Wx)):ceil(max(Wx)), floor(min(Wy)):ceil(max(Wy)));
i = inpolygon(denseWx, denseWy, Wx(boundary), Wy(boundary));
denseWx = denseWx(i); % remove points outside polygon
denseWy = denseWy(i);

% project coordinates into original image
denseOx = griddata(Wx, Wy, Ox, denseWx, denseWy, 'linear'); denseOx = denseOx(:);
denseOy = griddata(Wx, Wy, Oy, denseWx, denseWy, 'linear'); denseOy = denseOy(:);
%j=find(isnan(denseOx+denseOy)==0); 

% find points projected inside the target image
j = find(denseWx>0 & denseWx<ncols & denseWy>0 & denseWy<nrows);
denseOx=denseOx(j); denseWx=denseWx(j); denseOy=denseOy(j); denseWy=denseWy(j);


% image warping
[X,Y]=meshgrid(1:size(imgOriginal,2),1:size(imgOriginal,1));
for c=1:cc
    n = sub2ind([nrows, ncols, cc], denseWy, denseWx, c*ones(size(denseWx)));
    imgWarp(n) = uint8(interp2(X, Y, double(imgOriginal(:,:,c)), denseOx, denseOy,'linear'));
end
