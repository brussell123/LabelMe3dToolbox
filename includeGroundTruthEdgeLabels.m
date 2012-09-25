function annotation = includeGroundTruthEdgeLabels(annotation,edgeLabels,img)

CONTACT_EDGE = 1;
Nobjects = length(annotation.object);
largeObjects  = getListLargeObjects;
for i = 1:Nobjects
  if ismember(strtrim(lower(annotation.object(i).name)),largeObjects)
    Npoints = length(annotation.object(i).polygon.pt);
    edges = edgeLabels.object(i).edges;
    edgeType = edgeLabels.object(i).edgeType;
    n = find(edgeType==CONTACT_EDGE);
    c = edges(:,n);
    contactPoints = zeros(1,Npoints);
    contactPoints(unique(c(:))) = 1;
    for j = 1:Npoints
      annotation.object(i).polygon.pt(j).contact = contactPoints(j);
    end
    
% $$$     clf;
% $$$     imshow(img);
% $$$     hold on;
% $$$     [X,Y] = getLMpolygon(annotation.object(i).polygon);
% $$$     plot([X; X(1)],[Y; Y(1)],'b','LineWidth',4);
% $$$     for j = 1:length(X)
% $$$       if annotation.object(i).polygon.pt(j).contact
% $$$         plot(X(j),Y(j),'ro','LineWidth',4);
% $$$       end
% $$$     end
% $$$     pause;
  end
end
