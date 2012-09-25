function plotNmap(Nmap)
% Input:
% Nmap - uint32 map of object indices.

[nrows,ncols,cc] = size(Nmap);
maxVal = inf*ones(1,1,'uint32');
valid = find(Nmap(:)~=maxVal);
ndx = unique(Nmap(valid));

% Get colors for map:
colors = uint8(255*jet(double(max(ndx)+1)));

outR = zeros(nrows,ncols,'uint8');
outG = zeros(nrows,ncols,'uint8');
outB = zeros(nrows,ncols,'uint8');

for i = 1:length(ndx)
  nn = find(Nmap(valid)==ndx(i));
  outR(valid(nn)) = colors(ndx(i)+1,1);
  outG(valid(nn)) = colors(ndx(i)+1,2);
  outB(valid(nn)) = colors(ndx(i)+1,3);
end

outImg = reshape([outR outG outB],[nrows ncols 3]);
% $$$ colorbar('YTick',255*ytick,'YTickLabel',yticklabels,'YLim',255*ylim);
% $$$ colorbar('YTick',[0 double(max(ndx))],'YLim',[0 double(max(ndx))]);
% $$$ h = colorbar;
% $$$ set(h,'YLim',[0 double(max(ndx))]);

imshow(outImg);
