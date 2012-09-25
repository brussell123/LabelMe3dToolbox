% basic labels (all edges get one of these 5 labels):
% 1. contact edges
% 2. occlusion edge
% 3. occlusion edge, non-owner
% 4. attached edge (all non-occluded road/sidewalk edges are attached
% edges; we consider road as "part-of" ground plane)
% 5. other
%
% additional attributes:
%
% - visibility: all edges get one of these 3 labels:
% a. visible edge
% b. invisible edge
% c. partially invisible edge
%
% - for all contact edges, classify into one of these types:
% 1a. supported by ground
% 1b. supported by another object (include index pointing to object)
%
% - for all non-owner occlusion edges, indicate who owns the occlusion:
% 3a. image owns the occlusion (due to cropping)
% 3b. another object owns the occlusion (include index pointing to
% object)
%
% - for all attached edges, indicate who owns the part:
% 4a. which object (or ground plane) is this edge attached to (include
% index pointing to object)

% $$$ HOMEIMAGES = '/Users/brussell/work/Datasets/LabelMe/Images';
% $$$ LIBPATHS = {'/Users/brussell/work/MatlabLibraries/LabelMeToolbox'};
% $$$ 
% $$$ addpath(LIBPATHS{:});
% $$$ addpath('/Users/brussell/work/3Dsupport');
% $$$ addpath(genpath('/Users/brussell/work/3Dsupport/src'));
% $$$ 
% $$$ load /Users/brussell/work/3Dsupport/EdgeLabeling/DBbenchmark.mat;

% parameters:
resolutionTexture = 640;

labelingType = input('What type of labeling? [(b)asic, (v)isiblility] ','s');
switch labelingType
 case 'b'
  edgeOptions = {'Contact edge','Occlusion edge',sprintf('Occlusion edge\nnon-owner'),'Attached edge','Other edge','Start over'};
  labelingType = 'basic';
 case 'v'
  edgeOptions = {'Visible','partially-visible','invisible'};
  basicEdgeOptions = {'Contact edge','Occlusion edge',sprintf('Occlusion edge\nnon-owner'),'Attached edge','Other edge','Start over'};
  labelingType = 'visibility';
end

Njet = length(edgeOptions);
edgeTypeColors = jet;
edgeTypeColors = edgeTypeColors(round((size(edgeTypeColors,1)-1)*[0:Njet-1]/(Njet-1))+1,:);
for i = 13%:length(DB)
  display(sprintf('%d out of %d',i,length(DB)));
  
  annotation = DB(i).annotation;
  img = LMimread(DB,i,HOMEIMAGES);
  
  % re-scale and crop image
  [annotation,img] = LMimscale(annotation,img,min(1,resolutionTexture/size(img,1)));

  % Get number of edges:
  Nobjects = length(annotation.object);
  Nedges = 0;
  for j = 1:Nobjects
    [X,Y] = getLMpolygon(annotation.object(j).polygon);
    Nedges = Nedges+length(X);
  end

  j = 1;
  while j<=Nobjects
    [X,Y] = getLMpolygon(annotation.object(j).polygon);
    Nedges = Nedges-length(X);

    if exist('edgeLabels','var') && length(edgeLabels)>=i && length(edgeLabels(i).object)>=j && ((strcmp(labelingType,'basic') && isfield(edgeLabels(i).object(j),'edgeType') && ~any(edgeLabels(i).object(j).edgeType==0)) || (strcmp(labelingType,'visibility') && isfield(edgeLabels(i).object(j),'visibility') && ~isempty(edgeLabels(i).object(j).visibility) && ~any(edgeLabels(i).object(j).visibility==0)))
      % Display existing label:
      clf;
      subplot('position',[0 0 0.8 0.95]);
      plotEdgeLabels(annotation,img,edgeLabels(i),labelingType,j);
      title(sprintf('%s: %d edges, %d edges left',annotation.object(j).name,length(X),Nedges));
      
      % Query if user wants to edit edge labels:
      subplot('position',[0.85 0 0.15 1]);
      imagesc([1:Njet]');
      axis off;
      text(0.6,1,sprintf('Edit object label'),'Color',[1 1 1]);
      text(0.6,2,sprintf('Go to next label'),'Color',[1 1 1]);
      text(0.6,3,sprintf('Quit'));

      [x,y] = ginput(1);
      y = round(y)
      
      switch y
       case 2
        j = j+1;
        continue;
       case 3
        return;
      end
    end
    
    edges = [[1:length(X)]; [2:length(X) 1]];
    edgeType = zeros(1,length(X));
    
    k = 1;
    while k <= size(edges,2)
      if k==1
        clf;
        subplot('position',[0 0 0.8 0.95]);
        imshow(img);
        hold on;
        plot([X; X(1)],[Y; Y(1)],'w','LineWidth',4);
      end
      e = edges(:,k);
      subplot('position',[0 0 0.8 0.95]);
      plot([X(e(1)) X(e(2))],[Y(e(1)) Y(e(2))],'Color',edgeTypeColors(end,:),'LineWidth',4);
      switch labelingType
       case 'basic'
        title(sprintf('%s: %d edges, %d edges left',annotation.object(j).name,length(X),Nedges));
       case 'visibility'
        title(sprintf('%s (%s): %d edges, %d edges left',annotation.object(j).name,basicEdgeOptions{edgeLabels(i).object(j).edgeType(k)},length(X),Nedges));
      end

      subplot('position',[0.85 0 0.15 1]);
      imagesc([1:Njet]');
      axis off;

      for kk = 1:Njet
        text(0.6,kk,edgeOptions{kk},'Color',[1 1 1]);
      end
      
      [x,y] = ginput(1);
      y = round(y)

      if y==length(edgeOptions)
        k = 1;
        continue;
      end
      
      edgeType(k) = y;
      subplot('position',[0 0 0.8 0.95]);
      plot([X(e(1)) X(e(2))],[Y(e(1)) Y(e(2))],'Color',edgeTypeColors(y,:),'LineWidth',4);
      k = k+1;
    end
    
    % Plot all labels for object and confirm save:
    clf;
    subplot('position',[0 0 0.8 0.95]);
    imshow(img);
    hold on;
    for k = 1:size(edges,2)
      e = edges(:,k);
      plot([X(e(1)) X(e(2))],[Y(e(1)) Y(e(2))],'Color',edgeTypeColors(edgeType(k),:),'LineWidth',4);
    end

    subplot('position',[0.85 0 0.15 1]);
    imagesc([1:Njet]');
    axis off;
    text(0.6,1,sprintf('Save object label'),'Color',[1 1 1]);
    text(0.6,2,sprintf('Redo'),'Color',[1 1 1]);
    text(0.6,3,sprintf('Quit'));
    
    [x,y] = ginput(1);
    y = round(y)
    
    switch y
     case 1
      switch labelingType
       case 'basic'
        edgeLabels(i).object(j).edges = edges;
        edgeLabels(i).object(j).edgeType = edgeType;
       case 'visibility'
        edgeLabels(i).object(j).visibility = edgeType;
        if ~isfield(edgeLabels(i).object(j),'edges')
          edgeLabels(i).object(j).edges = edges;
        end
      end
      j = j+1;
     case 2
      Nedges = Nedges+length(X);
     case 3
      return;
    end
  end

  clf;
  plotEdgeLabels(annotation,img,edgeLabels(i),labelingType);

%  save edgeLabels.mat edgeLabels;
end

tot = 0;
for i = 1:length(edgeLabels)
  for j = 1:length(edgeLabels(i).object)
    tot = tot+size(edgeLabels(i).object(j).edges,2);
  end
end
tot

return;

% 6791

for i = 1:10
  annotation = DB(i).annotation;
  img = LMimread(DB,i,HOMEIMAGES);
  clf;
  plotEdgeLabels(annotation,img,edgeLabels(i),labelingType);
  print('-djpeg',sprintf('labeledImage%03d.jpg',i));
end
