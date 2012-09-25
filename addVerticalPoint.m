function [X,Y,CptsNdx] = addVerticalPoint(X,Y,xc,CptsNdx);
%
% Add points in the polygon that are the intersection of the edges of the
% polygon with the vertical lines defines by xc
% 
% Also, propagate the contact point label to the new points
%
% xc is a vector


bb = [min(X) max(X)];
delta = (bb(2)-bb(1))/100; 
nc = length(xc);
ndx = [];

for m = 1:nc
    X = [X; X(1)];
    Y = [Y; Y(1)];
    if nargin >3; CptsNdx = [CptsNdx; CptsNdx(1)]; end
    np = length(X)-1;
    % Search points (only add the point, if there is no point nearby)
    Xn = [];
    Yn = [];
    Cn = [];
    for i = 1:np
        if min(X(i),X(i+1))+delta<=xc(m) & max(X(i),X(i+1))-delta>xc(m)
            % Insert point
            xx = xc(m);
            yy = (Y(i+1)-Y(i))/(X(i+1)-X(i))*(xc(m)-X(i))+Y(i);
            Xn = [Xn; X(i); xx];
            Yn = [Yn; Y(i); yy];
            if nargin > 3
                Cn = [Cn; CptsNdx(i); CptsNdx(i)*CptsNdx(i+1)]; % the new point is a contact point if it is betwen two contact points
            end
        else
            Xn = [Xn; X(i)];
            Yn = [Yn; Y(i)];
            if nargin >3
                Cn = [Cn; CptsNdx(i)];
            end
        end
    end
    X = Xn;
    Y = Yn;
    if nargin >3
        CptsNdx = Cn;
    end
end
