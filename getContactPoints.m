function [x,y] = getContactPoints(polygon)
% [x,y] = getContactPoints(polygon)
%
% Get contact points for polygon.
%
% Inputs:
% polygon - LabelMe polygon structure.
%
% Outputs:
% x,y - Arrays of contact points.

x = [];
y = [];
if isfield(polygon,'contact')
  x = str2num(char({polygon.contact.x}));
  y = str2num(char({polygon.contact.y}));
end
