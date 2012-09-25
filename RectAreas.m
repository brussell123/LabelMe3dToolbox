function [int,union,A1,A2] = RectAreas(ro,rp)
% Fast intersection computation.  Used as a bound on the normalized
% intersection area.  Normalize by the area of the smaller rectangle.
% min_x,min_y,max_x,max_y  

  tt = 0;
  wo = ro(3)-ro(1);
  ho = ro(4)-ro(2);
  wp = rp(3)-rp(1);
  hp = rp(4)-rp(2);
  int = rectint([ro(1) ro(2) wo ho],[rp(1) rp(2) wp hp]);
  A1 = wo*ho;
  A2 = wp*hp;
  union = A1+A2-int;
