function [x,y] = samplePoly(X,Y,S)
% Evenly sample around polygon in clockwise fashion.

if nargin < 3
  S = 2000; % Number of sample points
end
x = zeros(S,1);
y = zeros(S,1);
N = length(X);

% Make input polygon clockwise:
if sum(X.*[Y(2:end); Y(1)]-[X(2:end); X(1)].*Y)<=0
  % Counter-clockwise, so flip order:
  X = flipud(X);
  Y = flipud(Y);
end
  
% Compute polygon spacing:
L = sum(sqrt((X-[X(2:end); X(1)]).^2+(Y-[Y(2:end); Y(1)]).^2))/S;

% Get top-mid point:
x(1) = (min(X)+max(X))/2;
c = find((X<=x(1))&(x(1)<=[X(2:end); X(1)]));
n = mod(c,N)+1;
a = (x(1)-X(n))./(X(c)-X(n));
[y(1),m] = max(a.*Y(c)+(1-a).*Y(n));

% Current, next, alpha, and distance vars:
c = c(m);
n = n(m);
a = 1-a(m);
d = sqrt((X(c)-X(n))^2+(Y(c)-Y(n))^2);

% Loop through polygon and get evenly spaced points:
for i = 2:S
  len = L;
  while len > 0
    l = a*d+len;
    if l<=d
      a = l/d;
      x(i) = (1-a)*X(c)+a*X(n);
      y(i) = (1-a)*Y(c)+a*Y(n);
      break;
    else
      len = l-d;
      c = n;
      n = mod(c,N)+1;
      a = 0;
      d = sqrt((X(c)-X(n))^2+(Y(c)-Y(n))^2);
    end
  end
end

return;

ff = figure;
plot([X; X(1)],[Y; Y(1)],'b','LineWidth',4);
hold on;
plot(x(1),y(1),'go','LineWidth',4);
plot(x(2),y(2),'ro','LineWidth',4);
plot(x(3:end),y(3:end),'bo','LineWidth',4);
axis ij;

pause;

close(ff);
