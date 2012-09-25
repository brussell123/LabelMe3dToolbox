function pb = isThereAnyPointBelow(x,y)
%
% A point can not be a contact point if there is another point in the same
% boundary that is bellow the current point and in the same vertical 
%
% this routine will return a indicator variable indicating if there is
% anything bellow it

x = [x;x(1)];
y = [y;y(1)];
np = length(x);
pb = zeros(np-1,1);

for n = 1:np-1
    yy = [];m=0;
    for i = 1:np-1
        if min(x(i),x(i+1))<=x(n)+0.001 & max(x(i),x(i+1))>x(n)-0.001 & i~=n
            % Insert point
            m = m+1;
            yy(m) = (y(i+1)-y(i))/(x(i+1)-x(i))*(x(n)-x(i))+y(i);
        end
    end
    if length(yy)>0
        pb(n) = y(n)<max(yy);
    end
end

