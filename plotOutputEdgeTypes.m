function plotOutputEdgeTypes(annotation,img,etype)

if nargin < 3
  etype = [1 2 3];
end

% $$$   edgeOptions = {'Attached edge','Contact edge','Occlusion edge'};
edgeOptions = {'Contact','Attached','Occlusion'};
Njet = length(edgeOptions);
edgeTypeColors = jet;
edgeTypeColors = edgeTypeColors(round((size(edgeTypeColors,1)-1)*[0:Njet-1]/(Njet-1))+1,:);

clf;
h = imshow(img);
set(h,'AlphaData',0.5);
hold on;
h = zeros(1,Njet);
for i = 1:length(annotation.object)
  [X,Y] = getLMpolygon(annotation.object(i).polygon);
  edges = annotation.object(i).polygon.edges;
  Nedges = size(edges,2);
  edgeType = annotation.object(i).polygon.edgeType;

  edgeType(~ismember(edgeType,etype)) = 0;
  edgeType = 2*(edgeType==1) + 1*(edgeType==2) + 3*(edgeType==3);
  for j = 1:Nedges
    if edgeType(j)>0
      e = edges(:,j);
      h(edgeType(j)) = plot([X(e(1)) X(e(2))],[Y(e(1)) Y(e(2))],'Color',edgeTypeColors(edgeType(j),:),'LineWidth',4);
    end
  end
end

if sum(h>0)>1
  legend(h(h>0),edgeOptions(h>0),'Location','NorthEastOutside');
else
  legend(edgeOptions(h>0),'Location','NorthEastOutside');
end
