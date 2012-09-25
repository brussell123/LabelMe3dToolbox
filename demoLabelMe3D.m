% This file contains the scripts needed to run the LabelMe3D algorithm.
% If you use this toolbox, please cite the following paper:
%
% B. C. Russell and A. Torralba. Building a Database of 3D Scenes from
% User Annotations. In CVPR, 2009.
%
% To start, change the paths below so that they point to the location of
% the LabelMe and LabelMe3D MATLAB Toolboxes.

% Add LabelMe and LabelMe3D MATLAB Toolboxes to your path:
addpath '~/work/MatlabLibraries/LabelMeToolbox';
addpath '~/work/MatlabLibraries/LabelMe3dToolbox';

% Set this if you want to write the outputs:
OUTDIR = '';

% LabelMe3D structures:
load Params3D.mat;

% Get image and LabelMe annotation:
fname_jpg = 'example.jpg';
fname_xml = 'example.xml';

img = imread(fname_jpg);
annotation = LMread(fname_xml);
  
% Recover 3D scene components:
imageSize = [size(img,1) size(img,2)];
[annotation,valid] = Recover3DSceneComponents(annotation,Params3D,imageSize);

if valid
  % Plot 3D scene:
  figure;
  LMplot3Dscene(annotation);

  % Show objects + horizon line:
  [xh,yh] = getHorizonLine(annotation);
  figure;
  LMplot(annotation,img);
  hold on;
  plot(xh,yh,'r');

  % Sanity check that 3D points project to 2D points:
  P = getCameraMatrix(annotation);
  err = zeros(2,length(annotation.object));
  for i = 1:length(annotation.object)
    if isfield(annotation.object(i).world3d,'polygon3d')
      X = getLMpolygon3D(annotation.object(i).world3d.polygon3d);
      x = P*[X; ones(1,size(X,2))];
      [x,y] = RH2LH(x(1,:)./x(3,:),x(2,:)./x(3,:),imageSize);
      [xx,yy] = getLMpolygon(annotation.object(i).polygon);
      err(1,:) = sum(abs(x'-xx));
      err(2,:) = sum(abs(y'-yy));
    end
  end
  display(sprintf('Mean error: %f',mean(err(:))));
  
  % Insert ground plane:
  annotation = insertGroundPlane(annotation);
  
% $$$   % Get XYZN maps:
% $$$   [Xmap,Ymap,Zmap,Nmap] = getXYZmaps(annotation);
% $$$   
% $$$   % Plot XYZ maps:
% $$$   figure;
% $$$   subplot(1,3,1); plotXmap(Xmap); title('X map');
% $$$   subplot(1,3,2); plotYmap(Ymap); title('Y map');
% $$$   subplot(1,3,3); plotZmap(Zmap); title('Z map');
% $$$   
% $$$   % Plot polygon and edge types:
% $$$   seg = plotPolyEdgeTypes(annotation,img,'valid');
% $$$   figure;
% $$$   imshow(seg);
% $$$   title('Polygon and edge types');

  % Get scene mesh:
  mesh = getSceneMesh(annotation);
  mesh = subDivideMesh(mesh,2);
  
  % Generate mesh textures:
  mesh = getTextures(mesh,annotation,img);

  figure;
  LMplot3Dmesh(mesh);
  
  if ~isempty(OUTDIR)
    % Write 3D XML file:
    writeXML3D(annotation,fullfile(OUTDIR,'example.out.xml'));
  
    % Create VRML file:
    vrmlfile = 'example.wrl';
    LM2VRML(annotation,mesh,vrmlfile,OUTDIR);
    
% $$$     % Create book
% $$$     BOOK = get3DBook(annotation,img);
% $$$     imwrite(BOOK,fullfile(OUTDIR,'example_book.jpg'),'jpg','quality',100);
% $$$ 
% $$$     % Write XYZN maps:
% $$$     writeXmap(Xmap,fullfile(OUTDIR,'example.X.png'));
% $$$     writeYmap(Ymap,fullfile(OUTDIR,'example.Y.png'));
% $$$     writeZmap(Zmap,fullfile(OUTDIR,'example.Z.png'));
% $$$     writeNmap(Nmap,fullfile(OUTDIR,'example.N.png'));
  end
end
