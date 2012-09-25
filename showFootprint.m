function showFootprint(annotation,img,HOME3D)

Hy = annotation.camera.Hy;
Nobjects = length(annotation.object);

% Plot ground objects on image:
clf;
imshow(img);
hold on;
plot([1 size(img,2)],[Hy Hy],'r','LineWidth',4);
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh) && strcmp(annotation.object(i).mesh.type,'ground')
    % Get coordinates:
    TX = annotation.object(i).mesh.v2D(1,:);
    TY = annotation.object(i).mesh.v2D(2,:);
    b = annotation.object(i).mesh.boundary;
    plot([TX(b) TX(b(1))],[TY(b) TY(b(1))],'b','LineWidth',4);
  end
end

if nargin > 2
  print('-djpeg',fullfile(HOME3D,annotation.folder,[annotation.filename '_groundObjects.jpg']));
else
  pause;
end

% Plot foldedcard objects on image:
clf;
imshow(img);
hold on;
plot([1 size(img,2)],[Hy Hy],'r','LineWidth',4);
colors = {'r','b','c','m'};
objColors = cell(1,Nobjects);
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh) && strcmp(annotation.object(i).mesh.type,'foldedcard') %&& sum(annotation.object(i).mesh.CptsNdx)>1
    TX = annotation.object(i).mesh.v2D(1,:);
    TY = annotation.object(i).mesh.v2D(2,:);
    b = annotation.object(i).mesh.boundary;
    rp = randperm(length(colors));
    objColors{i} = colors{rp(1)};
    plot([TX(b) TX(b(1))],[TY(b) TY(b(1))],objColors{i},'LineWidth',4);
  end
end
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh) && strcmp(annotation.object(i).mesh.type,'foldedcard') %&& sum(annotation.object(i).mesh.CptsNdx)>1
    TX = annotation.object(i).mesh.v2D(1,:);
    TY = annotation.object(i).mesh.v2D(2,:);
    b = annotation.object(i).mesh.boundary;
    c = find(annotation.object(i).mesh.CptsNdx);
    x = TX(b(c));
    y = TY(b(c));
    [v,j] = sort(x);
    if length(j) > 1
      plot(x(j),y(j),'g','LineWidth',4);
    else
      plot(x(j),y(j),'g.','MarkerSize',30);
    end      
  end
end

if nargin > 2
  print('-djpeg',fullfile(HOME3D,annotation.folder,[annotation.filename '_foldedCards.jpg']));
else
  pause;
end

% Plot footprint:
clf;
objName = [];
range = [inf -inf inf -inf];
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh)
    % Get coordinates:
    X = annotation.object(i).mesh.v3D(1,:);
    Y = annotation.object(i).mesh.v3D(2,:);
    Z = annotation.object(i).mesh.v3D(3,:);
    TX = annotation.object(i).mesh.v2D(1,:);
    TY = annotation.object(i).mesh.v2D(2,:);
    b = annotation.object(i).mesh.boundary;

    switch annotation.object(i).mesh.type
     case 'ground'
      plot([X(b) X(b(1))],[Z(b) Z(b(1))],'k','LineWidth',4);
      hold on;
      objName{end+1} = annotation.object(i).name;
      range([1 3]) = min(range([1 3]),[min(X) min(Z)]);
      range([2 4]) = max(range([2 4]),[max(X) max(Z)]);
     case 'foldedcard'
%      if sum(annotation.object(i).mesh.CptsNdx)>1
        x = X(b);
        z = Z(b);
% $$$         [x,j] = sort(x);
        [v,j] = sort(TX(b));
        x = x(j);
        z = z(j);
        plot(x,z,objColors{i},'LineWidth',4);
        hold on;
        objName{end+1} = annotation.object(i).name;
        range([1 3]) = min(range([1 3]),[min(X) min(Z)]);
        range([2 4]) = max(range([2 4]),[max(X) max(Z)]);
%      end
    end
  end
end
range(3) = min(range(3),0);
xlabel('x');
ylabel('z');
%axis equal;
axis(range,'equal');
grid on;
legend(objName,'Location','NorthEastOutside');

if nargin > 2
  print('-djpeg',fullfile(HOME3D,annotation.folder,[annotation.filename '_footprint.jpg']));
else
  pause;
end

return;

% Plot the remaining objects (foldedcards with a single contact point)
clf;
imshow(img);
hold on;
plot([1 size(img,2)],[Hy Hy],'r','LineWidth',4);
for i = 1:Nobjects
  if isfield(annotation.object(i),'mesh') && ~isempty(annotation.object(i).mesh) && strcmp(annotation.object(i).mesh.type,'foldedcard') && sum(annotation.object(i).mesh.CptsNdx)==1
    TX = annotation.object(i).mesh.v2D(1,:);
    TY = annotation.object(i).mesh.v2D(2,:);
    b = annotation.object(i).mesh.boundary;
    plot([TX(b) TX(b(1))],[TY(b) TY(b(1))],'b','LineWidth',4);
  end
end

if nargin > 2
  print('-djpeg',fullfile(HOME3D,annotation.folder,[annotation.filename '_singleContactPoint.jpg']));
else
  pause;
end
