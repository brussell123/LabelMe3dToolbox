function yc = intersectPolygon(x,y,xc)
% finds all the vertical intersection points 'yc' 
% of the polygon [x,y]

yc = [];
x = [x; x(1)];
y = [y; y(1)];

n = length(x);

j = find((x(1:n-1)<=xc & x(2:n)>=xc) |  (x(2:n)<=xc & x(1:n-1)>=xc) );
for k = 1:length(j)
    i = j(k);
    if x(i+1)==x(i) 
        yc = [yc y(i+1) y(i)]; % vertically aligned, so add both points
    else
        yc = [yc (y(i+1)-y(i))/(x(i+1)-x(i))*(xc-x(i))+y(i)];
    end
end
