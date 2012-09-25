function [pd,pf,areaROC,Scores] = roc(Scores,GT)

GT = GT>0;

[Scores,i] = sort(Scores,'descend');
GT = GT(i);

pd = cumsum(GT);
pf = cumsum(~GT);

% Get unique scores:
[s,i] = unique(Scores,'last');
i = i(end:-1:1);
pd = pd(i)/sum(GT);
pf = pf(i)/sum(~GT);
Scores = Scores(i);

pd = [0 pd 1];
pf = [0 pf 1];
Scores = [inf Scores -inf];

k = unique(convhull([pf 1],[pd 0]));
pd = pd(k(1:end-1));
pf = pf(k(1:end-1));
Scores = Scores(k(1:end-1));

areaROC = polyarea([pf 1],[pd 0]);
