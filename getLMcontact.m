function [x,y] = getLMcontact(polygon)
% Inputs:
% polygon
% 
% Outputs:
% x
% y

x = []; y = [];
if isfield(polygon,'contact')
  x = str2num(char({polygon.contact.x}));
  y = str2num(char({polygon.contact.y}));
end
