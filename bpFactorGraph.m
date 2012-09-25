function beliefs = bpFactorGraph(nodes,factors)
% Run belief propagation given factor graph.
%
% Inputs:
% nodes(i).msgIn - Messages from factors.
% nodes(i).msgOut - Messages to factors.
%
% factors(i).pot;
% factors(i).ndxNodes;
% factors(i).ndxMsgs;

% Parameters:
Niterations = 10;

Nnodes = length(nodes);
Nfactors = length(factors);
for i = 1:Niterations
  % Pass messages from nodes to factors:
  for j = 1:Nnodes
    msg = log(nodes(j).msgIn);
    m = msg;
    m(isinf(msg)) = 0;
    msg = repmat(sum(msg,2),1,size(msg,2))-m;

    % Normalize messages:
    msg = exp(msg-repmat(max(msg,[],1),size(msg,1),1));
    msg = msg./repmat(sum(msg,1),size(msg,1),1);

    % Set outgoing messages:
    nodes(j).msgOut = msg;
  end
  
  % Pass messages from factors to nodes:
  for j = 1:Nfactors
    pot = factors(j).pot;
    ndxNodes = factors(j).ndxNodes;
    ndxMsgs = factors(j).ndxMsgs;
    
    switch ndims(pot)
     case 2
      if prod(size(pot))==length(pot)
        nodes(ndxNodes).msgIn(:,ndxMsgs) = pot;
      else
        m1 = nodes(ndxNodes(1)).msgOut(:,ndxMsgs(1));
        m2 = nodes(ndxNodes(2)).msgOut(:,ndxMsgs(2));
        m = pot*m2;
        nodes(ndxNodes(1)).msgIn(:,ndxMsgs(1)) = m/sum(m);
        m = m1'*pot;
        nodes(ndxNodes(2)).msgIn(:,ndxMsgs(2)) = m/sum(m);
      end
     case 3
      m1 = nodes(ndxNodes(1)).msgOut(:,ndxMsgs(1));
      m2 = nodes(ndxNodes(2)).msgOut(:,ndxMsgs(2));
      m3 = nodes(ndxNodes(3)).msgOut(:,ndxMsgs(3));
      m = zeros(size(pot,1),1);
      for s = 1:size(pot,1)
        m(s) = m2'*squeeze(pot(s,:,:))*m3;
      end
      nodes(ndxNodes(1)).msgIn(:,ndxMsgs(1)) = m/sum(m);
      m = zeros(size(pot,2),1);
      for s = 1:size(pot,2)
        m(s) = m1'*squeeze(pot(:,s,:))*m3;
      end
      nodes(ndxNodes(2)).msgIn(:,ndxMsgs(2)) = m/sum(m);
      m = zeros(size(pot,3),1);
      for s = 1:size(pot,3)
        m(s) = m1'*squeeze(pot(:,:,s))*m2;
      end
      nodes(ndxNodes(3)).msgIn(:,ndxMsgs(3)) = m/sum(m);
     otherwise
      error('BP can only handle 2 or 3 node factors for now.');
    end
  end
end

% Compute beliefs:
for j = 1:Nnodes
  msg = sum(log(nodes(j).msgIn),2);
  
  % Normalize messages:
  msg = exp(msg-max(msg));
  beliefs(j).bel = msg/sum(msg);
end
