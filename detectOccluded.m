function occluded = detectOccluded(annotation,img,Ppart,objectClasses)
% Detect which objects are occluded.

Nobjects = length(annotation.object);
[nrows,ncols,dim] = size(img);
  
% First, detect parts so that they don't occlude the object they belong
% to:
parts = inferParts(annotation,Ppart,objectClasses);

% Sort layers:
[a,i,layers] = LMsortlayers(annotation,img);
% $$$ a == annotation(i)

% Get polygon masks:
for j = 1:Nobjects
  [X{j},Y{j}] = getLMpolygon(annotation.object(j).polygon);
% $$$   [X{j},Y{j}] = samplePoly(X{j},Y{j});
end

occluded = logical(zeros(1,Nobjects));
for j = 1:Nobjects
  for k = 1:Nobjects
    if (j~=k) && (layers(i==k)>=layers(i==j)) && (parts(k)~=j) && ...
           (any(inpolygon(X{j},Y{j},X{k},Y{k})) || ...
            any(inpolygon(X{k},Y{k},X{j},Y{j})))
      occluded(j) = 1;
      break;
    end
  end
end

return;

figure;
for j = 1:length(occluded)
  if ~occluded(j) && ~parts(j)
    clf;
    imshow(img);
    hold on;
    plot([X{j}; X{j}(1)],[Y{j}; Y{j}(1)],'b','LineWidth',4);
    pause;
  end
end

im_counts = LMcountobject(DB);
[v,n] = sort(im_counts,'descend');
figure;
c = 1;
imnum = 1105;%1079;%1052;%1035;%1018;%1000;
while c<=64
  annotation = DB(n(imnum)).annotation;
  img = LMimread(DB,n(imnum),HOMEIMAGES);
  sup = inferSupport(annotation,Psup,objectClasses);
  unocc = find(~detectOccluded(annotation,img,Ppart,objectClasses) & ...
               ~inferParts(annotation,Ppart,objectClasses) & ...
               sup);
  unocc = unocc(sup(sup(unocc))==0);
  
  j = 1;
  while (j<=length(unocc))&&(c<=64)
    imgCrop = LMobjectcrop(img,annotation,unocc(j),2);
    subplottight(8,8,c);
    imshow(imgCrop);
    drawnow;
    j = j+1;
    c = c+1;
  end
  imnum = imnum+1
end
