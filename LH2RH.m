function [x,y] = LH2RH(x,y,imageSize)

x = (imageSize(2)+1)/2-x;
y = (imageSize(1)+1)/2-y;
% $$$ x = x-(imageSize(2)+1)/2;
% $$$ y = (imageSize(1)+1)/2-y;

