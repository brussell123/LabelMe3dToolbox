function Hy = findHorizonLine(annotation,img)

% For now, have the user select the horizon line.  This will need to be
% found automatically.
figure;
imshow(img);
title('Select horizon line');
[Hx Hy] = ginput(1);  % get horizon from user
