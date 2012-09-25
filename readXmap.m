function Xmap = readXmap(fname)
  
% Read from PNG file

png = int32(imread(fname));
sgn = int32(-2*floor(double(png(:,:,3))/128)+1);
Xmap = zeros(size(png,1),size(png,2),'int32');
Xmap = png(:,:,1)+256*png(:,:,2)+256*256*mod(png(:,:,3),128);
Xmap = Xmap.*sgn;
