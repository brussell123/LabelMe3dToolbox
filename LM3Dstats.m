function [countObj,countImg] = LM3Dstats(HOMEANNOTATIONS)

DB = LM3Ddatabase(HOMEANNOTATIONS);

countImg = 0;
countObj = 0;
for i = 1:length(DB)
  annotation = DB(i).annotation;
  cc = 0;
  for j = 1:length(annotation.object)
    if isfield(annotation.object(j),'world3d') && isfield(annotation.object(j).world3d,'polygon3d')
      countObj = countObj+1;
      cc = 1;
    end
  end
  countImg = countImg+cc;
end
