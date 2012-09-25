function [x,y] = RH2LH(x,y,imageSize)

x = (imageSize(2)+1)/2-x;
y = (imageSize(1)+1)/2-y;
% $$$ x = x+(imageSize(2)+1)/2;
% $$$ y = (imageSize(1)+1)/2-y;
