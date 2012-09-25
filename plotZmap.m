function Zmap = plotZmap(Zmap,lim,units)
% Zmap = plotZmap(Zmap,lim)
%
% Output pretty depth map figure.  If no output argument is specified,
% then a figure is created showing the depth map
%
% Inputs:
% Zmap - Depth map in absolute coordinates.
% lim - Min and max depth limits (optional)
% units - Use 'metric' or 'english' units (optional, 'metric' is default)
%
% Outputs:
% Zmap - Output depth map, with colors indicating absolute distance
  
MAXVAL = 16777215;
if (nargin < 2) || isempty(lim)
  %  MAX_DEPTH = 16777215; % max uint24
  MIN_DEPTH = 1e2; % 1m
  MAX_DEPTH = 1e5; % 1km
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
Zinvalid = find(Zmap==MAXVAL);

% Make Zmap lie within MIN_DEPTH/MAX_DEPTH range:
Zmap = log10(min(max(double(Zmap),MIN_DEPTH),MAX_DEPTH));

% Adjust Zmap range to lie between 0 and 1:
Zmap = Zmap-ylim(1);
ytick = ytick-ylim(1);
ylim = ylim-ylim(1);
Zmap = Zmap/ylim(2);
ytick = ytick/ylim(2);
ylim = ylim/ylim(2);

% Pass Zmap through -jet colormap:
cm = uint8(255*fliplr(jet(64)));
z = round(Zmap*(size(cm,1)-1))+1;
Zmap(:,:,1) = cm(z);
Zmap(:,:,2) = cm(z+size(cm,1));
Zmap(:,:,3) = cm(z+2*size(cm,1));
Zmap = uint8(Zmap);

% Make invalid regions black:
Zmap(Zinvalid) = 0;
Zmap(Zinvalid+size(Zmap,1)*size(Zmap,2)) = 0;
Zmap(Zinvalid+2*size(Zmap,1)*size(Zmap,2)) = 0;

if ~nargout
%  figure;
  imagesc(Zmap,[0 255]);
%  imagesc(Zmap,[0 1]);
  colormap(double(cm)/255);
  axis image off;
  colorbar('YTick',255*ytick,'YTickLabel',yticklabels,'YLim',255*ylim);
end
