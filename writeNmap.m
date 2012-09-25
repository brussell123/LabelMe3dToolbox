function writeNmap(Nmap,fname)
  
% Write to PNG file

MAX_DEPTH = 16777215;

Nmap = max(min(Nmap,MAX_DEPTH),0);

png = zeros(size(Nmap,1),size(Nmap,2),3,'uint8');
png(:,:,1) = mod(Nmap,256);
png(:,:,2) = mod(bitshift(Nmap,-8),256);
png(:,:,3) = mod(bitshift(Nmap,-16),256);

imwrite(png,fname,'png');
