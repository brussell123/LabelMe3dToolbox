function [Xmap,Ymap,Zmap,Nmap] = getXYZmaps(annotation,mesh,scaleFactor)

if nargin < 3
% $$$   maxRes = 360;
  scaleFactor = 1;
end
if nargin < 2
  mesh = getSceneMesh(annotation);
  mesh = subDivideMesh(mesh,2);
end

P = getCameraMatrix(annotation);

imageSize = [str2num(annotation.imagesize.nrows) str2num(annotation.imagesize.ncols)];

if scaleFactor~=1%max(imageSize)>maxRes
% $$$   scaleFactor = maxRes/max(imageSize);
  P = [scaleFactor 0 0; 0 scaleFactor 0; 0 0 1]*P;
  imageSize = round(scaleFactor*imageSize);
end

[K,R,C] = decomposeP(P);

% Ray-cast:
[x,y] = meshgrid(1:imageSize(2),1:imageSize(1));
[x,y] = LH2RH(x,y,imageSize);
x = [x(:) y(:) ones(prod(size(x)),1)]';
[X,isValid,faceNdx] = rayTrace(single(P(:,1:3)\x),single(C),mesh.vertices,mesh.faces);

isValid = reshape(isValid,imageSize);
faceNdx = reshape(faceNdx,imageSize);
Nmap = inf*ones(imageSize,'uint32');
Nmap(isValid) = mesh.objndx(faceNdx(isValid));
Xmap = int32(reshape(X(1,:),imageSize));
Ymap = uint32(reshape(X(2,:),imageSize));
Zmap = uint32(reshape(X(3,:),imageSize));

MAX_DEPTH = 16777215;
MAX_WIDTH = 8388607;
Xmap(~isValid) = MAX_WIDTH;
Ymap(~isValid) = MAX_DEPTH;
Zmap(~isValid) = MAX_DEPTH;

% $$$ function [Xmap,Ymap,Zmap,Nmap,Xobj,Yobj,Zobj] = getXYZmaps_old(annotation)
% $$$ % [Xmap,Ymap,Zmap,Nmap] = getXYZmaps(annotation)
% $$$ %
% $$$ % Computes depth map given the computed 3D mesh.
% $$$ % Returns the X,Y,Z data for each pixel in the image, along with a
% $$$ % map of object IDs (corresponding to annotation.object(i).id)
% $$$ % indicating which object is closest to the camera center at each x,y
% $$$ % location in the image.  Pixels that do not point to any valid object
% $$$ % with 3D information is given the value of 4294967295.
% $$$ %
% $$$ % [Xmap,Ymap,Zmap,Nmap,Xobj,Yobj,Zobj] = getXYZmaps(annotation)
% $$$ %
% $$$ % Computes individual depth maps for the objects.
% $$$ %
% $$$ % Inputs:
% $$$ %   annotation - LabelMe annotation structure
% $$$ %
% $$$ % Outputs:
% $$$ %   Xmap - Map of X coordinates
% $$$ %   Ymap - Map of Y coordinates
% $$$ %   Zmap - Map of Z coordinates
% $$$ %   Nmap - Map of object IDs for each pixel (ground and standing objects).
% $$$ %   Xobj - 3D matrix of X coordinates for each object.
% $$$ %   Yobj - 3D matrix of Y coordinates for each object.
% $$$ %   Zobj - 3D matrix of Z coordinates for each object.
% $$$ 
% $$$ MAX_DEPTH = 16777215;
% $$$ MAX_WIDTH = 8388607;
% $$$ nrows = str2num(annotation.imagesize.nrows);
% $$$ ncols = str2num(annotation.imagesize.ncols);
% $$$ 
% $$$ % $$$ % Insert ground plane:
% $$$ % $$$ annotation = insertGroundPlane(annotation);
% $$$ 
% $$$ % Get indices of non-deleted objects:
% $$$ notDeleted = find(~isdeleted(annotation))';
% $$$ 
% $$$ Nobjects = length(annotation.object);
% $$$ Xmap = MAX_WIDTH*ones(nrows,ncols,'int32');
% $$$ Ymap = MAX_DEPTH*ones(nrows,ncols,'uint32');
% $$$ Zmap = MAX_DEPTH*ones(nrows,ncols,'uint32');
% $$$ Nmap = inf*ones(nrows,ncols,'uint32');
% $$$ if nargout>=5
% $$$   Xobj = MAX_WIDTH*ones(nrows,ncols,Nobjects,'int32');
% $$$   Yobj = MAX_DEPTH*ones(nrows,ncols,Nobjects,'uint32');
% $$$   Zobj = MAX_DEPTH*ones(nrows,ncols,Nobjects,'uint32');
% $$$ end
% $$$ 
% $$$ % Store depth maps for standing and ground objects:
% $$$ for i = notDeleted%1:Nobjects
% $$$   if ismember(annotation.object(i).polygon.polyType,{'ground','standing'})
% $$$     [X,Y,Z,mx,my,tri,boundary,CptsNdx,valid] = getObject3D(annotation,i,'grid');
% $$$     
% $$$     if ~isempty(X) && valid
% $$$       X = int32(X)';
% $$$       Y = uint32(Y)';
% $$$       Z = uint32(Z)';
% $$$       n = sub2ind([nrows ncols],my,mx);
% $$$       nn = ((Z<Zmap(n)) | (Nmap(n)==inf*ones(1,1,'uint32')));
% $$$       Xmap(n(nn)) = X(nn);
% $$$       Ymap(n(nn)) = Y(nn);
% $$$       Zmap(n(nn)) = Z(nn);
% $$$       Nmap(n(nn)) = str2num(annotation.object(i).id);
% $$$       if nargout>=5
% $$$         Xobj(n+(i-1)*nrows*ncols) = X;
% $$$         Yobj(n+(i-1)*nrows*ncols) = Y;
% $$$         Zobj(n+(i-1)*nrows*ncols) = Z;
% $$$       end
% $$$     end
% $$$   end
% $$$ end
% $$$ 
% $$$ % Store depth maps for attached objects:
% $$$ for i = notDeleted%1:Nobjects
% $$$   if strcmp(annotation.object(i).polygon.polyType,'part')
% $$$     [X,Y,Z,mx,my,tri,boundary,CptsNdx,valid] = getObject3D(annotation,i,'grid');
% $$$     parents = getPartOfParents(annotation,i);
% $$$ 
% $$$     if ~isempty(X) && valid && ~isempty(parents)
% $$$       X = int32(X)';
% $$$       Y = uint32(Y)';
% $$$       Z = uint32(Z)';
% $$$       n = sub2ind([nrows ncols],my,mx);
% $$$       nn = ((Z<Zmap(n)) | (Nmap(n)==inf*ones(1,1,'uint32')));
% $$$       Xmap(n(nn)) = X(nn);
% $$$       Ymap(n(nn)) = Y(nn);
% $$$       Zmap(n(nn)) = Z(nn);
% $$$       Nmap(n(nn)) = str2num(annotation.object(i).id);
% $$$       nn = [];
% $$$       for j = 1:length(parents)
% $$$         nn = [nn find(Nmap(n)==str2num(annotation.object(parents(j)).id))];
% $$$       end
% $$$       Nmap(n(nn)) = str2num(annotation.object(i).id);
% $$$       if nargout>=5
% $$$         Xobj(n+(i-1)*nrows*ncols) = X;
% $$$         Yobj(n+(i-1)*nrows*ncols) = Y;
% $$$         Zobj(n+(i-1)*nrows*ncols) = Z;
% $$$       end
% $$$     end
% $$$   end
% $$$ end
% $$$ 
% $$$ if nargout==0
% $$$   plotDepthMap(Zmap);
% $$$ end
