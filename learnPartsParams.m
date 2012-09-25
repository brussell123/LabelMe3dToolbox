function [par,N] = learnPartsParams(evidence,objectClasses)
% Inputs:
% evidence - Structure of evidence.
% objectClasses - Cell array of object classes.
%
% Outputs:
% P(i,j) - Probability that object class i is attached to object class
% j.  If j=length(objectClasses)+1, then this is probability of not being
% attached to any object class.

% Parameters:
Niterations = 10; % Number of Gibbs sampling iterations
alpha_r = 20; % Beta distribution prior for relative overlap
alpha_a = 2;%0.8; % Beta distribution prior for relative area
eta_a = 5; % Beta prior for no attachment
eta_o = 0.5; % Dirichlet prior for which object attached to

% Get number of objects and seed random number generator:
Nobjects = length(objectClasses);
rand('twister',sum(100*clock));

% Collect data point information:
imgNdx = [];
objNdx = [];
for i = 1:length(evidence)
  if ~isempty(evidence(i).objNdx)
    n = find(evidence(i).objNdx>0);
    imgNdx = [imgNdx i*ones(1,length(n),'int32')];
    objNdx = [objNdx int32(n)];
  end
end

% Run Gibbs sampler:
Ndata = length(imgNdx);
N = zeros(Nobjects,Nobjects+2,'int32'); % Last two cols are detect counts
q = -1*ones(1,Ndata,'int32');
for i = 1:Niterations
  display(sprintf('%d out of %d',i,Niterations));
  rp = int32(randperm(Ndata));
  [N,q] = partsGibbsSampler(evidence,imgNdx,objNdx,N,q,rp,alpha_r,alpha_a,eta_a,eta_o);
end

[Ndata sum(sum(N(:,1:end-1))) sum(sum(N(:,end-1:end)))]

% Compute output parameters:
P = double(N);
P(:,end-1) = (P(:,end-1)+eta_a)./(P(:,end-1)+P(:,end)+eta_a);
P(:,1:end-2) = P(:,1:end-2)+eta_o;
P(:,1:end-2) = P(:,1:end-2)./repmat(sum(P(:,1:end-2),2),1,Nobjects);
P(:,end) = [];

par.P = P;
par.alpha_r = alpha_r;
par.alpha_a = alpha_a;
par.eta_a = eta_a;
par.eta_o = eta_o;

return;

%%%% TEST %%%%

[v,n] = sort(Ppart.P(:,end));
for i = 1:200
  display(sprintf('%s: %f',objectClasses{n(i)},1-v(i)));
  [v2,n2] = sort(Ppart.P(n(i),1:Nobjects),'descend');
  for j = 1:5
    if v2(j) > 10*eta_o/Nobjects
      display(sprintf('-> %s: %f',objectClasses{n2(j)},v2(j)));
    end
  end
  pause;
end

% license plate: 435
[v,n] = sort(P(435,1:Nobjects),'descend');

% Visualize beta distribution
figure;
alpha_a = 10;%0.8;
x = [0:0.01:1];
p = alpha_a*x.^(alpha_a-1);
plot(x,p);
axis([0 1 0 max(p)]);

% Collect evidence statistics:
c = 793;
for i = 1:length(evidence)
  objNdx = evidence(i).objNdx;
  relativeArea = evidence(i).relativeArea;
  n = find(objNdx==c);
  for j = 1:length(n)
    relativeOverlap = evidence(i).relativeOverlap(n(j),:);
    [v,k] = sort(relativeOverlap,'descend');
    v(k==n(j)) = [];
    k(k==n(j)) = [];
    display(sprintf('%s: %f; %s: %f',objectClasses{objNdx(k(1))},v(1),objectClasses{objNdx(k(2))},v(2)));
    pause;
  end
end
