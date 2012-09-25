function plotEdgeLabels(annotation,img,edgeLabels,labelingType,i)

Nobjects = length(annotation.object);

if nargin < 5
  i = 1:Nobjects;
end

switch labelingType
 case 'basic'
  edgeOptions = {'Contact edge','Occlusion edge',sprintf('Occlusion edge\nnon-owner'),'Attached edge','Other edge'};
 case 'visibility'
  edgeOptions = {'Visible','partially-visible','invisible'};
 case 'simplified'
  edgeOptions = {'Contact','Attached','Occlusion'};
end

Njet = length(edgeOptions);
edgeTypeColors = jet;
edgeTypeColors = edgeTypeColors(round((size(edgeTypeColors,1)-1)*[0:Njet-1]/(Njet-1))+1,:);

imshow(img);
hold on;
h = zeros(1,Njet);
for j = i
  edges = edgeLabels.object(j).edges;
  switch labelingType
   case 'basic'
    edgeType = edgeLabels.object(j).edgeType;
   case 'visibility'
    edgeType = edgeLabels.object(j).visibility;
   case 'simplified'
    edgeType = edgeLabels.object(j).edgeType;
    edgeType(ismember(edgeType,[2 3 5])) = 3;
    edgeType(edgeType==4) = 2;
  end
  [X,Y] = getLMpolygon(annotation.object(j).polygon);
  for k = 1:size(edges,2)
    e = edges(:,k);
    h(edgeType(k)) = plot([X(e(1)) X(e(2))],[Y(e(1)) Y(e(2))],'Color',edgeTypeColors(edgeType(k),:),'LineWidth',4);
  end
end

if sum(h>0)>1
  legend(h(h>0),edgeOptions(h>0),'Location','NorthEastOutside');
else
  legend(edgeOptions(h>0),'Location','NorthEastOutside');
end
