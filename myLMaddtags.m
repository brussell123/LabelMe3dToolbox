function annotation = myLMaddtags(annotation,tagsfile)

if nargin < 2
  tagsfile = 'tags.txt';
end

% Below is code extracted from LMaddtags.m

D.annotation = annotation;

[Tag, Descriptions] = loadtags(tagsfile);
Ntags = length(Tag);

[labelmeDescriptions,counts,imagendx,objectndx,descriptionndx] = LMobjectnames(D);
ndxtag = zeros(length(labelmeDescriptions),1);
for i = 1:length(labelmeDescriptions)
  for k = 1:Ntags
    j = strmatch(lower(labelmeDescriptions{i}), Descriptions{k}, 'exact');
    if ~isempty(j)
      ndxtag(i) = k;
      break
    end
  end
end
for i = 1:length(labelmeDescriptions)
  if ndxtag(i)>0
    j = find(descriptionndx==i);
    for k = 1:length(j)
      D(imagendx(j(k))).annotation.object(objectndx(j(k))).name = Tag{ndxtag(i)};
    end
  else
    j = find(descriptionndx==i);
    for k = 1:length(j)
      D(imagendx(j(k))).annotation.object(objectndx(j(k))).name = 'unmatched';
    end
  end
end

for i = 1:length(annotation.object)
  annotation.object(i).originalname = annotation.object(i).name;
  annotation.object(i).name = D.annotation.object(i).name;
end
