function Ymap = readYmap(fname)
  
% Read from PNG file

png = uint32(imread(fname));
Ymap = zeros(size(png,1),size(png,2),'uint32');
Ymap = png(:,:,1)+256*png(:,:,2)+256*256*png(:,:,3);
