function valid = isCleanObject(X,Y,imageSize)
% valid = isCleanObject(X,Y,imageSize)
%
% Returns true if the polygon is not cropped and is not very small.
%
% Inputs:
% X - X coordinates of polygon
% Y - Y coordinates of polygon
% imageSize - [nrows ncols] of image
%
% Outputs:
% valid - 1 if polygon satisfies conditions.

% Parameters:
minHeightPixels = 0.02;%0.1; % Minimum height of object in pixels
cropDist = 0.01; % Minimum distance to image boundary before considered
                 % as cropped object

minHeightPixels = round(imageSize(1)*minHeightPixels);
cropDist = (imageSize(1)^2+imageSize(2)^2)^0.5*cropDist;

valid = 0;
if ((max(Y)-min(Y))>=minHeightPixels) && (min(Y)>cropDist) && (min(X)>cropDist) && (max(Y)<(imageSize(1)-cropDist)) && (max(X)<(imageSize(2)-cropDist))
  valid = 1;
end
