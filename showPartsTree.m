function showPartsTree(Ppart, partObjectClasses, object, DB)
%
% showPartsTree(Ppart, partObjectClasses)
%
% showPartsTree(Ppart, partObjectClasses, object)

th = .75;

ndx = strmatch(lower(object), lower(partObjectClasses), 'exact');

% search children (parts of object)
children = find(Ppart(ndx,:)>th);
[foo,j] = sort(Ppart(ndx,children), 'descend');
children = children(j);

% search (object belong to)
parents = find(Ppart(:,ndx)>th);
[foo,j] = sort(Ppart(parents,ndx), 'descend');
parents = parents(j);

N = max(length(parents), length(children));
xo = 1; yo = N/2;

[foo,n] = sort(Ppart(ndx,:), 'descend');
for k=children
    disp(sprintf('%s is part of %s (with probability = %1.2f)', partObjectClasses{k}, object, Ppart(ndx,k)))
end

figure
hold on
axis([-1 3 0 N+1])
for c = 1:length(children)
    xc = 2; yc = N/2-(length(children)+1)/2+c;
    plot([xo xc], [yo yc], 'k')
end
for c = 1:length(parents)
    xc = 0; yc = N/2-(length(parents)+1)/2+c;
    plot([xo xc], [yo yc], 'k')
end
for c = 1:length(children)
    xc = 2; yc = N/2-(length(children)+1)/2+c;
    text (xc, yc, strrep(partObjectClasses{children(c)}, '_', ' '), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', [1 1 0], 'EdgeColor','k')
end
for c = 1:length(parents)
    xc = 0; yc = N/2-(length(parents)+1)/2+c;
    text (xc, yc, strrep(partObjectClasses{parents(c)}, '_', ' '), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', [1 .2 0], 'EdgeColor','k')
end
text (xo, yo, object, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', [1 .7 0], 'EdgeColor','k')
axis('off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 3
    % Show examples
    [foo,j] = LMquery(DB, 'object.name', object);
    [foo,jp] = LMquery(DB(j), 'object.name', partObjectClasses{children(1)});

    figure
    for m = 1:min(4, length(jp))
        annotation = DB(j(jp(m))).annotation;

        % Get object mask
        img = LMimread(DB, j(jp(m)), 'http://labelme.csail.mit.edu/Images');
        [nrows, ncols, cc] = size(img);
        [annotation,img] = LMimscale(annotation, img, 320/ncols);
        [nrows, ncols] = getaproximagesize(annotation);
        mask = LMobjectmask(annotation, [ncols nrows], object);

        subplot(2,4,2*m-1)
        LMplot(DB, img);
        legend off
        
        subplot(2,4,2*m)
        imagesc(1-mask); colormap(gray(256)); axis('equal'); axis('off')
        hold on

        % Get part polygons
        colors = hsv(length(children));
        for p = 1:length(children)
            jc = LMobjectindex(annotation, partObjectClasses{children(p)});
            for k = 1:length(jc)
                [X,Y] = getLMpolygon(annotation.object(jc(k)).polygon);
                plot([X; X(1)], [Y; Y(1)], 'LineWidth', 2, 'color', colors(p,:))
            end
        end
    end
end




