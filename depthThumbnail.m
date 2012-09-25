function thumbnail = depthThumbnail(annotation, HOMEIMAGES, img, Zmap)
%
% To call the function, first inflate a scene:
%  [annotation,Zmap]=InflateScene(DB(ndx).annotation,Hy,sup,parts,Cpts,img);
%  thumbnail=depthThumbnail(DB(ndx).annotation, HOMEIMAGES, [], Zmap);
%  figure
%  imshow(thumbnail) 


MIN_DEPTH = 2e2; % 1m
MAX_DEPTH = 2e4; % 1km

thumb = LMsceneThumbnail(annotation, HOMEIMAGES, img, 'lines');
[nrows ncols cc] = size(thumb)

thumb = thumb(:, fix(ncols/2)+2:ncols,:);


mask = Zmap<500000;
Zmap = log10(min(max(double(Zmap),MIN_DEPTH),MAX_DEPTH));
Zmap = uint8(256*(Zmap-log10(MIN_DEPTH)) / (log10(MAX_DEPTH) - log10(MIN_DEPTH)));

Zmap = Zmap.*uint8(mask);
cmap = fliplr(jet(256));
cmap(1,:) = [0 0 0];

depthmap = ind2rgb(Zmap, cmap);
depthmap = imresize(depthmap, [size(thumb,1) size(thumb,2)], 'bilinear');
depthmap = uint8(255*depthmap);

thumbnail = [thumb 255*ones(nrows,3,3) depthmap];
