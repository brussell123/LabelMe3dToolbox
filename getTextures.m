function [mesh,tn] = getTextures(mesh,annotation,img)
% Gets textures for the objects.
%
% Inputs:
% mesh - Triangular mesh (note: mesh.vertices has 3xN entries)
% annotation - LabelMe annotation structure
% img - Image
%
% Outputs:
% mesh.textures - Cell array of object textures.
% mesh.tx - 2xN texture coordinates corresponding to 3D vertices.
% tn - 1xN texture indices
%
%
% Example:
%
% mesh = getTextures(mesh,annotation,img);
% 
% Display mesh points on texture:
%
% i - Desired object index.
% n = find(mesh.objndx==i);
% figure;
% imshow(mesh.textures{i});
% hold on;
% triplot(mesh.faces(:,n)',mesh.tx(1,:),mesh.tx(2,:));

P = getCameraMatrix(annotation);
[K,R,C] = decomposeP(P);
imageSize = [size(img,1) size(img,2)];

% Ray-cast:
[x,y] = meshgrid(1:imageSize(2),1:imageSize(1));
[x,y] = LH2RH(x,y,imageSize);
x = [x(:) y(:) ones(prod(size(x)),1)]';
[X,isValid,faceNdx] = rayTrace(single(P(:,1:3)\x),single(C),mesh.vertices,mesh.faces);

% Get object map:
isValid = reshape(isValid,imageSize);
faceNdx = reshape(faceNdx,imageSize);
objMap = zeros(imageSize);
objMap(isValid) = mesh.objndx(faceNdx(isValid));

% Get textures:
N = size(mesh.vertices,2);
objNdx = unique(mesh.objndx);
mesh.textures = cell(1,length(annotation.object));
mesh.tx = zeros(2,N,'single');
tn = zeros(1,N);
for i = objNdx
  % Project mesh points:
  n = mesh.faces(:,find(mesh.objndx==i));
  n = double(unique(n(:)))';
  xi = P*[mesh.vertices(:,n); ones(1,length(n))];
  [xi,yi] = RH2LH(xi(1,:)./xi(3,:),xi(2,:)./xi(3,:),imageSize);
  
  % Get bounding box:
  bb = [max(1,floor(min(xi))) min(imageSize(2),ceil(max(xi))) max(1,floor(min(yi))) min(imageSize(1),ceil(max(yi)))];

  % Get crops:
  imgCrop = img(bb(3):bb(4),bb(1):bb(2),:);
  objMapCrop = objMap(bb(3):bb(4),bb(1):bb(2));
  
% $$$   figure;
% $$$   imshow(imgCrop);
% $$$   figure;
% $$$   imagesc(objMapCrop);
  
  % Get average color within object:
  imgCrop_r = double(squeeze(imgCrop(:,:,1)));
  imgCrop_g = double(squeeze(imgCrop(:,:,2)));
  imgCrop_b = double(squeeze(imgCrop(:,:,3)));
  colorAvg = uint8([mean(imgCrop_r(objMapCrop==i)) mean(imgCrop_g(objMapCrop==i)) mean(imgCrop_b(objMapCrop==i))]);

  % Fill non-object pixels with average color:
  imgCrop_r(objMapCrop~=i) = colorAvg(1);
  imgCrop_g(objMapCrop~=i) = colorAvg(2);
  imgCrop_b(objMapCrop~=i) = colorAvg(3);
  imgCrop = uint8(reshape([imgCrop_r imgCrop_g imgCrop_b],size(imgCrop)));

  % Get texture coordinates:
  mesh.tx(:,n) = [xi-bb(1); yi-bb(3)];
  tn(n) = i;
  mesh.textures{i} = imgCrop;
  
% $$$   figure;
% $$$   imshow(imgCrop);
% $$$   hold on;
% $$$   plot(xi,yi,'r.');
end
