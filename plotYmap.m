function Ymap = plotYmap(Ymap,lim,units)
% Ymap = plotDepthMap(Ymap,lim)
%
% Output pretty height map figure.  If no output argument is specified,
% then a figure is created showing the height map
%
% Inputs:
% Ymap - Height map in absolute coordinates.
% lim - Min and max height limits (optional)
% units - Use 'metric' or 'english' units (optional, 'metric' is default)
%
% Outputs:
% Ymap - Output height map, with colors indicating absolute distance
  
MAXVAL = 16777215;
if (nargin < 2) || isempty(lim)
  %  MAX_DEPTH = 16777215; % max uint24
  MIN_DEPTH = 1e1;%1e2; % 1m
  MAX_DEPTH = 1e4;%1e5; % 1km
else
  MIN_DEPTH = lim(1);
  MAX_DEPTH = lim(2);
end
if nargin < 3
  units = 'metric';
end

ytick = [0:7];
switch units
 case 'metric'
  yticklabels = {'1cm','10cm','1m','10m','100m','1km','10km','100km'};
 case 'english'
  yticklabels = {'0.03ft','0.33ft','3.28ft','32.81ft','328.08ft','3280.80ft','32808ft','100km'};
 otherwise
  error('Invalid units');
end

ylim = [log10(MIN_DEPTH) log10(MAX_DEPTH)];

% Get valid regions:
Zinvalid = find(Ymap==MAXVAL);

% Make Ymap lie within MIN_DEPTH/MAX_DEPTH range:
Ymap = log10(min(max(double(Ymap),MIN_DEPTH),MAX_DEPTH));

% Adjust Ymap range to lie between 0 and 1:
Ymap = Ymap-ylim(1);
ytick = ytick-ylim(1);
ylim = ylim-ylim(1);
Ymap = Ymap/ylim(2);
ytick = ytick/ylim(2);
ylim = ylim/ylim(2);

% Pass Ymap through -jet colormap:
cm = uint8(255*fliplr(jet(64)));
z = round(Ymap*(size(cm,1)-1))+1;
Ymap(:,:,1) = cm(z);
Ymap(:,:,2) = cm(z+size(cm,1));
Ymap(:,:,3) = cm(z+2*size(cm,1));
Ymap = uint8(Ymap);

% Make invalid regions black:
Ymap(Zinvalid) = 0;
Ymap(Zinvalid+size(Ymap,1)*size(Ymap,2)) = 0;
Ymap(Zinvalid+2*size(Ymap,1)*size(Ymap,2)) = 0;

if ~nargout
%  figure;
  imagesc(Ymap,[0 255]);
%  imagesc(Ymap,[0 1]);
  colormap(double(cm)/255);
  axis image off;
  colorbar('YTick',255*ytick,'YTickLabel',yticklabels,'YLim',255*ylim);
end
