function [z F] = avg_error(Z, G);
% Z is our guess
% G is ground truth

% Parameters:
F = 0;

[nrows,ncols,d] = size(Z);

% rescale G (note that all images in this dataset have the same size):
GT = imresize(G, [nrows ncols]);


E = 0;
count = 0;
for i = 1:nrows
  for j = 1:ncols
    Z(i,j) = double(Z(i,j))/100;
    if (Z(i,j) > 0)
      % just to avoid divide-by-zeros
      count = count+1;
      if (GT(i,j) == 0)
        GT(i,j) = 0.1;
      end;
      E = E + abs(double(Z(i,j))-GT(i,j))/GT(i,j);
      F = F + double(Z(i,j))-GT(i,j);
    end;
  end;
end;
z = E/count;

return;

% $$$   % Z is our guess
% $$$   % G is ground truth
% $$$   F = 0;
% $$$ 
% $$$   s = size(Z);
% $$$   h = s(1);
% $$$   w = s(2);
% $$$ 
% $$$   %rescale G, since it has a funky dimension
% $$$   % note that all images in this dataset have the same size
% $$$   GT = imresize(G, s);
% $$$   
% $$$ 
% $$$   E = 0;
% $$$   count = 0;
% $$$   for i = 1:h
% $$$     for j = 1:w
% $$$       Z(i,j) = double(Z(i,j))/100;
% $$$       if (Z(i,j) > 0)
% $$$ 
% $$$ 	% just to avoid divide-by-zeros
% $$$ 	count = count+1;
% $$$ 	if (GT(i,j) == 0)
% $$$ 	  GT(i,j) = 0.1;
% $$$ 	end;
% $$$ 	E = E + abs(double(Z(i,j))-GT(i,j))/GT(i,j);
% $$$ 	F = F + double(Z(i,j))-GT(i,j);
% $$$ 	%	E = E + abs(double(Z(i,j))-GT(i,j))/GT(i,j);
% $$$       end;
% $$$     end;
% $$$   end;
% $$$   z = E/count;
