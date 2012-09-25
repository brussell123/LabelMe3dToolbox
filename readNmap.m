function Nmap = readNmap(fname)
  
% Read from PNG file

png = uint32(imread(fname));
Nmap = zeros(size(png,1),size(png,2),'uint32');
Nmap = png(:,:,1)+256*png(:,:,2)+256*256*png(:,:,3);
