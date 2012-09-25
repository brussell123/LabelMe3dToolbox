function annotation = inferEdgeTypesNew(annotation,validObjects)
% annotation = inferEdgeTypesNew(annotation)
%
% Infer edge types for all the labeled polygons in the scene.  Edges can
% either be (i) attached, (ii) occluded, (iii) contact.
%
% Inputs:
% annotation
% validObjects - Vector of desired object indices to compute.
%
% Outputs:
% annotation.object(i).polygon.edgeType
% annotation.object(i).polygon.pEdgeType
%
% Edge types:
% attached - 1
% contact  - 2
% occluded - 3

% Parameters:
priorContact = 0.05; % Prior of being a contact (or occluded) edge
beta_g = 50; % Beta distribution prior for distance to ground object
beta_l = 1.5; % Beta distribution prior for edge length
beta_o = 2; % Beta distribution prior for edge orientation
rho1 = 0.01; % Trinary term weight
rho2 = 10; % Trinary term weight
minEdgeLen = 0.05; % Minimum contact edge length
insertEdgeTypes = 0;

if nargin < 2
  validObjects = 1:length(annotation.object);
end

% Collect edge evidence:
[evidenceEdges,annotation] = collectEdgeEvidence(annotation,validObjects);

Nobjects = length(annotation.object);
ncols = annotation.imagesize.ncols;
for i = validObjects
  if strcmp(annotation.object(i).polygon.polyType,'standing')
    edges = annotation.object(i).polygon.edges;
    Nedges = size(edges,2);
    pEdgeType = zeros(3,Nedges);
    edgeType = zeros(1,Nedges);
    
    % Get edge evidence:
    pBelow = evidenceEdges.object(i).pBelow;
    pOri = evidenceEdges.object(i).pOri;
    pLen = evidenceEdges.object(i).pLen;
    pGround = evidenceEdges.object(i).pGround;
    
    % Get attached log-probability:
    pAttached = annotation.object(i).polygon.pAttached;
    pEdgeType(1,:) = log(pAttached);
    
    % Look for additional edge constraints:
    concaveConstraints = detectShallowConcave(annotation.object(i).polygon,pBelow);
    
    % Set up factor graph:
    % 1 - contact; 2 - occlusion
    k = 0;
    clear nodes factors;
    for j = 1:Nedges
      % Local edge evidence:
      k = k+1;
      p = [priorContact*betapdf(pGround(j),beta_g,1)*betapdf(pOri(j),beta_o,1)*betapdf(pLen(j),beta_l,1)*pBelow(j) 1-priorContact]';
      factors(k).pot = p/sum(p);
      factors(k).ndxNodes = j;
      factors(k).ndxMsgs = 1;
      nodes(j).msgIn = p;
    end
    
    if ~isempty(concaveConstraints)
      for j = 1:size(concaveConstraints,1)
        k = k+1;
        c = concaveConstraints(j,:);
        p = ones(2,2,2);
        p(2,1,2) = rho1; % contact
        p(2,2,2) = rho2; % occlusion
        factors(k).pot = p;
        factors(k).ndxNodes = c;
        factors(k).ndxMsgs = [size(nodes(c(1)).msgIn,2)+1 size(nodes(c(2)).msgIn,2)+1 size(nodes(c(3)).msgIn,2)+1];
        nodes(c(1)).msgIn(:,end+1) = ones(2,1);
        nodes(c(2)).msgIn(:,end+1) = ones(2,1);
        nodes(c(3)).msgIn(:,end+1) = ones(2,1);
      end
    end
    
    % Run BP:
    beliefs = bpFactorGraph(nodes,factors);
    
    % Assign beliefs:
    pEdgeType = [beliefs(:).bel];
    
    % Assign edge type based on maximal probability:
    [v,edgeType] = max(pEdgeType,[],1);
    
    % Set contact edges:
    edges = annotation.object(i).polygon.edges;
    Npts = length(annotation.object(i).polygon.pt);
    [X,Y] = getLMpolygon(annotation.object(i).polygon);
    
    % Get contact edges:
    c = zeros(1,Npts);
    n = find(edgeType==1);
    e = edges(:,n);
    c(unique(e(:))) = 1;

    % If only one edge segment and it is short, then use horizontal card:
    if sum(c)==2
      n = find(c);
      edgeLen = ((X(n(1))-X(n(2))).^2+(Y(n(1))-Y(n(2))).^2).^0.5/ncols;
      if edgeLen < minEdgeLen
        [v,j] = min([Y(n(1)) Y(n(2))]);
        c(n(j)) = 0;
      end
    end
    
    % Get contact points:
    n = find(c);
    annotation.object(i).polygon = setLMcontact(X(n),Y(n),annotation.object(i).polygon);
  end

  % Insert edge types (optional):
  if insertEdgeTypes
    for j = 1:length(annotation.object(i).polygon.pt)
      if strcmp(annotation.object(i).polygon.polyType,'standing')
        switch edgeType(j)
         case 1
          annotation.object(i).polygon.pt(j).edgeType = 'c';
         case 2
          annotation.object(i).polygon.pt(j).edgeType = 'o';
        end
      else
        annotation.object(i).polygon.pt(j).edgeType = 'a';
      end
    end
  end
  
end

return;

figure;
x = [0:0.01:1];
plot(x,betapdf(x,1.5,1));
axis([0 1 0 1]);

for i = 1:length(annotation.object)
  [X,Y] = getLMpolygon(annotation.object(i).polygon);
  edges = annotation.object(i).polygon.edges;  
  concaveConstraints = detectShallowConcave(annotation.object(i).polygon);
  if ~isempty(concaveConstraints)
    for j = 1:size(concaveConstraints,1)
      clf;
      imshow(img);
      hold on;
      plot([X; X(1)],[Y; Y(1)],'b','LineWidth',4);
      c = concaveConstraints(j,:);
      plot([X(edges(1,c(1))) X(edges(2,c(1)))],[Y(edges(1,c(1))) Y(edges(2,c(1)))],'r','LineWidth',4);
      plot([X(edges(1,c(2))) X(edges(2,c(2)))],[Y(edges(1,c(2))) Y(edges(2,c(2)))],'r','LineWidth',4);
      plot([X(edges(1,c(3))) X(edges(2,c(3)))],[Y(edges(1,c(3))) Y(edges(2,c(3)))],'r','LineWidth',4);
      pause;
    end
  end
end
