function LM2VRMLfill(annotation,vrmlfolder,outfile,imfile,img)
% LM2VRML(annotation,outfile,imfile,img)
%
% Converts LabelMe 3D polygons into VRML file for visualization.
%
% Inputs:
%   annotation - LabelMe annotation structure
%
% outfile - Output VRML file name (e.g. 'foo.wrl').
%   imfile - Pointer to location of image filename.  This is used to
%            texturize 3D polygons.
%   img -    Image.  If this is an empty matrix, then the VRML output will not
%            have image texture (the polygons will be colored).

% Get image size
if nargin>=4
    [nrows,ncols,dim] = size(img);
else
    if ~isfield(annotation,'imagesize')
        error('Please pass the image as an input or set the imagesize field.');
        return;
    end
    nrows = annotation.imagesize.nrows;
    ncols = annotation.imagesize.ncols;
    img = 1;
end

% Parameters VRML:
fieldOfView = 0.78;
scaling = 30;
vBB = [-10^5 10^5 -10^5 10^5 -100 10^8]; % volume bounding box (remove objects outside this volume)

% Get meshes
[X,Y,Z,TX,TY,tri,boundary,textures,objNames] = getTextures(annotation, img);
Nobjects = length(X);

% Save textures 
relativeFolder = strrep(outfile, '.wrl', '');
mkdir(fullfile(vrmlfolder, relativeFolder))
for ii = 1:Nobjects
    texturefilename{ii} = sprintf('object_%03d.jpg', ii);
    imwrite(textures{ii}, fullfile(vrmlfolder, relativeFolder, texturefilename{ii}), 'jpg')
end

% Get initial camera position:
initialCamera = getInitialViewpoint(X,Y,Z);

% Dump VRML file:
fp = fopen(fullfile(vrmlfolder,outfile),'w');
fprintf(fp,'#VRML V2.0 utf8\n\nNavigationInfo {\n  headlight TRUE\n  type [\"EXAMINE\", \"ANY\"]\n}\n\n');
fprintf(fp,'Viewpoint {\n  position %.2f %.2f %.2f\n  orientation 0 0 1 0\n  fieldOfView %.2f\n  description "Original"\n}\n\n',initialCamera(1), initialCamera(2), initialCamera(3), fieldOfView);
fprintf(fp,'DEF Back1 Background {\n  groundColor [1 1 1]\n  skyColor [1 1 1]\n}\n\n');

for ii = 1:Nobjects
    fprintf(fp,'# %s\nShape{\n',objNames{ii});
    if isempty(img)
        fprintf(fp,'  appearance Appearance {\n    material Material {\n      diffuseColor %.2f %.2f %.2f\n    }\n  }\n',colors(ii,1),colors(ii,2),colors(ii,3));
    else
        fprintf(fp,'  appearance Appearance {\n    texture ImageTexture { url \"%s\" }\n  }\n', fullfile(relativeFolder, texturefilename{ii}));
    end

    fprintf(fp,'  geometry IndexedFaceSet {\n    coord Coordinate {\n      point [\n');
    for jj = 1:length(X{ii})
        fprintf(fp,'        %.4f %.4f %.4f,\n',X{ii}(jj)/scaling,Y{ii}(jj)/scaling,-Z{ii}(jj)/scaling);
    end
    fprintf(fp,'      ]\n    }\n');

    %% Get coordIndex:
    fprintf(fp,'    coordIndex [\n');
    for j = 1:size(tri{ii},1)
      fprintf(fp,'      %d %d %d -1,\n',tri{ii}(j,:));
    end
    fprintf(fp,'    ]\n');

    if ~isempty(img)
        %% Get texture coordinates:
        fprintf(fp,'    texCoord TextureCoordinate {\n      point [\n');
        for jj = 1:length(TX{ii})
            fprintf(fp,'        %.4f %.4f,\n',TX{ii}(jj),1-TY{ii}(jj));
        end
        fprintf(fp,'      ]\n');
        fprintf(fp,'    }\n');

        %% Get texture coordinate indices:
        fprintf(fp,'    texCoordIndex [\n');
        for j = 1:size(tri{ii},1)
          fprintf(fp,'      %d %d %d -1,\n',tri{ii}(j,:));
        end
        fprintf(fp,'    ]\n');
    end

    fprintf(fp,'    solid FALSE\n');
    fprintf(fp,'  }\n');
    fprintf(fp,'}\n\n');
end
fclose(fp);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialCamera = getInitialViewpoint(X,Y,Z)
%
% From the object coordinater, get best view point to start viewing the
% scene

% Get initial camera position:
cube = [Inf -Inf Inf -Inf Inf -Inf 0 0 0];
for ii = 1:length(X)
    %cube = [minX maxX minY maxY minZ maxZ];
    cube(1) = min(cube(1), min(X{ii}));
    cube(2) = max(cube(2), max(X{ii}));
    cube(3) = min(cube(3), min(Y{ii}));
    cube(4) = max(cube(4), max(Y{ii}));
    cube(5) = min(cube(5), min(Z{ii}));
    cube(6) = max(cube(6), max(Z{ii}));
    % mean 
    cube(7) = cube(7)+mean(X{ii})/length(X);
    cube(8) = cube(8)+mean(max(1,Y{ii}))/length(X);
    cube(9) = cube(9)+mean(Z{ii})/length(X); 
end
%initialCamera = [(cube(1)+cube(2))/2  max(1,(cube(3)+cube(4))/2) cube(5)*0.8];
initialCamera = [cube(7) cube(8) 0*cube(5)*0.8];

initialCamera = [0 5 max(-30,min(cube(5),0))];


