function seg = plotLineBmap(seg,edges,color,sz)
% Plots lines on a bitmap image.

if nargin < 4
  sz = 2;
end

N = ceil(max(((edges(:,3)-edges(:,1)).^2+(edges(:,4)-edges(:,2)).^2).^0.5))+10;
alpha = linspace(0,1,N);
sx = zeros(size(edges,1),N);
sy = zeros(size(edges,1),N);
for i = 1:N
  sx(:,i) = (1-alpha(i))*edges(:,1)+alpha(i)*edges(:,3);
  sy(:,i) = (1-alpha(i))*edges(:,2)+alpha(i)*edges(:,4);
end
s = round([sx(:) sy(:)]);
s = max(s,1);
s(:,1) = min(s(:,1),size(seg,2));
s(:,2) = min(s(:,2),size(seg,1));
s = unique(s,'rows');
ndx = sub2ind(size(seg),s(:,2),s(:,1));
mask = zeros(size(seg,1),size(seg,2));
mask(ndx) = 1;
cB = strel('disk',sz);
mask = imdilate(mask,cB);
ndx = find(mask(:));

seg(ndx) = round(255*color(1));
seg(ndx+size(seg,1)*size(seg,2)) = round(255*color(2));
seg(ndx+2*size(seg,1)*size(seg,2)) = round(255*color(3));

