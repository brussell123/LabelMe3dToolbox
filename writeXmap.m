function png = writeXmap(Xmap,fname)
  
% Write to PNG file

MAX_WIDTH = 8388607;

Xmap = max(min(Xmap,MAX_WIDTH),-MAX_WIDTH);
Xmap = double(Xmap);
XmapSign = 1-min(sign(Xmap)+1,1);
Xmap = abs(Xmap);

png = zeros(size(Xmap,1),size(Xmap,2),3,'uint8');
png(:,:,1) = mod(Xmap,256);
png(:,:,2) = mod(bitshift(Xmap,-8),256);
png(:,:,3) = mod(bitshift(Xmap,-16),256)+128*XmapSign;

imwrite(png,fname,'png');
