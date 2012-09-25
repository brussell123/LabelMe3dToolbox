function [evidence,annotation] = collectEdgeEvidence(annotation,validObjects)
% Inputs:
% annotation
% validObjects - Vector of desired object indices to compute.
%
% Outputs:
% annotation.object(i).polygon.edges
% evidence.object(i).pBelow
% evidence.object(i).pOri
% evidence.object(i).pLen
% evidence.object(i).pGround

% Parameters:
Nsamples = 10;

if nargin < 2
  validObjects = 1:length(annotation.object);
end

Nobjects = length(annotation.object);
if isfield(annotation,'imagesize')
  nrows = annotation.imagesize.nrows;
  ncols = annotation.imagesize.ncols;
else
  [ncols,nrows] = getaproximagesize(annotation);
end

notDeleted = find(~isdeleted(annotation))';

% Get ground object indices:
groundObjects = getListGroundObjects; 
groundObjects = sprintf('%s,', groundObjects{:}); groundObjects = groundObjects(1:end-1);
groundObjects = LMobjectindex(annotation, groundObjects, 'exact');
groundObjects = groundObjects(ismember(groundObjects,notDeleted));

% Get polygons and add edges:
for i = notDeleted
  [X{i},Y{i}] = getLMpolygon(annotation.object(i).polygon);
end

evidence = [];
alpha = linspace(1,0,Nsamples);
for i = validObjects
  if strcmp(annotation.object(i).polygon.polyType,'standing')
    if ~isfield(annotation.object(i).polygon,'edges')
      annotation.object(i).polygon.edges = [[1:length(X{i})]; [2:length(X{i}) 1]];
    end

    otherGroundObjects = setdiff(groundObjects,i);
    edges = annotation.object(i).polygon.edges;
    Nedges = size(edges,2);
    
    % Do not consider as contact edge if there is a point below:
    pBelow = 1-double(isThereAnyPointBelowEdge(X{i},Y{i},edges));
    
    % Measure edge verticality:
    pOri = 1-edgeVerticality(X{i},Y{i},edges)';
    
    % Measure edge length width respect to object width:
    len = edgeLength(X{i},Y{i},edges)';
    objWidth = max(X{i})-min(X{i});
    pLen = min(len/objWidth,1);
    
    % Compute distance to ground object:
    pGround = zeros(1,Nedges);
    for j = 1:Nedges
      % Evenly sample over edge:
      e = edges(:,j);
      x = alpha*X{i}(e(1))+(1-alpha)*X{i}(e(2));
      y = alpha*Y{i}(e(1))+(1-alpha)*Y{i}(e(2));
      dGround = inf*ones(1,Nsamples);
      
      % Compute distance of edge samples to edges on other polygons:
      for ii = notDeleted
        if any(ii==otherGroundObjects)
          d = zeros(1,Nsamples);
          for jj = 1:Nsamples
            d(jj) = min(SupportDist(X{ii},Y{ii},x(jj),y(jj))/ncols,1);
          end
          dGround = min(d,dGround);
        end
      end
      
      % Compute probability of edge type by taking mean over edge
      % samples:
      pGround(j) = 1-mean(dGround);
    end
    
    % Store evidence:
    evidence.object(i).pBelow = pBelow;
    evidence.object(i).pOri = pOri;
    evidence.object(i).pLen = pLen;
    evidence.object(i).pGround = pGround;
  end
end
