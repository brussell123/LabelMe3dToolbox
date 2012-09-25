function Xmap = plotXmap(Xmap,lim,units)
% Xmap = plotXmap(Xmap,lim)
%
% Output pretty width map figure.  If no output argument is specified,
% then a figure is created showing the width map
%
% Inputs:
% Xmap - Width map in absolute coordinates.
% lim - Min and max width limits (optional)
% units - Use 'metric' or 'english' units (optional, 'metric' is default)
%
% Outputs:
% Xmap - Output width map, with colors indicating absolute distance

MAXVAL = 8388607;%16777215;
if (nargin < 2) || isempty(lim)
  v = double(max(abs(Xmap(Xmap~=MAXVAL))));
  v = 1000;
  MIN_DEPTH = -v;
  MAX_DEPTH = v;
else
  MIN_DEPTH = lim(1);
  MAX_DEPTH = lim(2);
end
if nargin < 3
  units = 'metric';
end

ytick = round(linspace(MIN_DEPTH,MAX_DEPTH,7));
yticklabels = [];
for i = 1:length(ytick)
  switch units
   case 'metric'
    yticklabels{i} = sprintf('%dcm',ytick(i));
   case 'english'
    yticklabels{i} = sprintf('%0.2fft',ytick(i)*unitsratio('feet','cm'));
   otherwise
    error('Invalid units');
  end
end

% Get valid regions:
Zinvalid = find(Xmap==MAXVAL);

Xmap = min(max(double(Xmap),MIN_DEPTH),MAX_DEPTH);

% Adjust Xmap range to lie between 0 and 1:
Xmap = Xmap-MIN_DEPTH;
ytick = ytick-MIN_DEPTH;
Xmap = Xmap/(MAX_DEPTH-MIN_DEPTH);
ytick = ytick/(MAX_DEPTH-MIN_DEPTH);

% Pass Xmap through -jet colormap:
cm = uint8(255*fliplr(jet(64)));
z = round(Xmap*(size(cm,1)-1))+1;
Xmap(:,:,1) = cm(z);
Xmap(:,:,2) = cm(z+size(cm,1));
Xmap(:,:,3) = cm(z+2*size(cm,1));
Xmap = uint8(Xmap);

% Make invalid regions black:
Xmap(Zinvalid) = 0;
Xmap(Zinvalid+size(Xmap,1)*size(Xmap,2)) = 0;
Xmap(Zinvalid+2*size(Xmap,1)*size(Xmap,2)) = 0;

if ~nargout
  imagesc(Xmap,[0 255]);
  colormap(double(cm)/255);
  axis image off;
  colorbar('YTick',255*ytick,'YTickLabel',yticklabels);
end
