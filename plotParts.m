function plotParts(annotation,img,parts)

% Visualize support relations:
Nobjects = length(annotation.object);
parents = find(parts==0);
Nsubplot = length(unique(parts))+2;
Nx = ceil(sqrt(Nsubplot));
Ny = ceil(Nsubplot/Nx);

%% Resize image and annotation for space savings:
[annotation,img] = LMimscale(annotation,img,min(640/size(img,2),1));
for i = 1:Nobjects
  [X{i},Y{i}] = getLMpolygon(annotation.object(i).polygon);
end

figure;
subplot(Ny,Nx,1);
LMplot(annotation,img);
legend hide;

k = 2;
i = 1;
while i <= length(parents)
  if ismember(parents(i),parts)
    subplot(Ny,Nx,k);
    imshow(img);
    hold on;
    Xi = X{parents(i)};
    Yi = Y{parents(i)};
    plot([Xi; Xi(1)],[Yi; Yi(1)],'r','LineWidth',4);
    children = find(parts==parents(i));
    for j = 1:length(children)
      Xi = X{children(j)};
      Yi = Y{children(j)};
      plot([Xi; Xi(1)],[Yi; Yi(1)],'b','LineWidth',4);
    end
    parents = [parents children(~ismember(children,parents))];
    title('Part relations');
    k = k+1;
  end
  i = i+1;
end

% $$$ subplot(Ny,Nx,k);
% $$$ showTreeObj(parts,{annotation.object(:).name});
% $$$ axis off;
% $$$ title('Support tree');
