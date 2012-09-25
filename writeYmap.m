function writeYmap(Ymap,fname)
  
% Write to PNG file

MAX_DEPTH = 16777215;

Ymap = max(min(Ymap,MAX_DEPTH),0);

png = zeros(size(Ymap,1),size(Ymap,2),3,'uint8');
png(:,:,1) = mod(Ymap,256);
png(:,:,2) = mod(bitshift(Ymap,-8),256);
png(:,:,3) = mod(bitshift(Ymap,-16),256);

imwrite(png,fname,'png');
