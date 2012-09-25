function polygon = setLMcontact(x,y,polygon)
% Inputs:
% x
% y
% polygon
%
% Outputs:
% polygon

% Need to have unique x, else interp1 does not work.  Need to handle
% this in a better way.
[x,j] = unique(x);
y = y(j);

for i = 1:length(x)
  polygon.contact(i).x = num2str(x(i));
  polygon.contact(i).y = num2str(y(i));
end
