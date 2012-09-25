function [R1,R2] = relativeOverlap(X1,Y1,X2,Y2)
% r = relativeOverlap(X1,Y1,X2,Y2)
% Compute relative overlap between two polygons.
%
% Inputs:
% X1,Y1 - Polygon 1.
% X2,Y2 - Polygon 2.
%
% Outputs:
% R1 - Relative overlap of polygon 1: area(P1 \cup P2)/area(P1).
% R2 - Relative overlap of polygon 2: area(P1 \cup P2)/area(P2).

% Do fast approximation of intersection area:
Rn = [min(X1) min(Y1) max(X1) max(Y1)];
Rm = [min(X2) min(Y2) max(X2) max(Y2)];
if Rn(3)<Rm(1) || Rm(3)<Rn(1) || Rn(4)<Rm(2) || Rm(4)<Rn(2)
  R1 = 0;
  R2 = 0;
  return;
end

% Compute relative overlap:
[A,ua,area1,area2] = PolyAreas(X1,Y1,X2,Y2);
R1 = A/(area1+eps);
R2 = A/(area2+eps);
