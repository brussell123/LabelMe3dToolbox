function writeZmap(Zmap,fname)
  
% Write to PNG file

MAX_DEPTH = 16777215;

Zmap = max(min(Zmap,MAX_DEPTH),0);

png = zeros(size(Zmap,1),size(Zmap,2),3,'uint8');
png(:,:,1) = mod(Zmap,256);
png(:,:,2) = mod(bitshift(Zmap,-8),256);
png(:,:,3) = mod(bitshift(Zmap,-16),256);

imwrite(png,fname,'png');
