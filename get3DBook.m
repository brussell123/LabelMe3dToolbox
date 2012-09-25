function BOOK = get3DBook(annotation, img)
% 
% Create all the information needed to build the VRML file

% TO DO:
%   - constraint the volume of the scene and crop the book. Maybe
%   contraint the farthest point


booksize = 1024*2;
maxZ = 1000; % remove objects farther than this threshold.


% Get meshes
[X,Y,Z,TX,TY,tri,boundary,textures,objNames,annotation] = getTextures(annotation, img);
Nobjects = length(X);

% Reduce size of 3D coordinates:
for n = 1:Nobjects
  X{n} = X{n}/10;
  Y{n} = Y{n}/10;
  Z{n} = Z{n}/10;
end

% Volume bounding box
vbb = [Inf -Inf Inf -Inf]; % [minZ maxZ minX maxX]
for n = 1:Nobjects
    if strcmp(annotation.object(n).polygon.polyType, 'ground')
        j = find(Z{n}<=maxZ);

        if ~isempty(j)
            vbb(1) = min(vbb(1), min(Z{n}(j)));
            vbb(2) = max(vbb(2), max(Z{n}(j)));
            vbb(3) = min(vbb(3), min(X{n}(j)));
            vbb(4) = max(vbb(4), max(X{n}(j)));
        end
    else
        CptsNdx = boundary{n}(find(annotation.object(n).mesh.CptsNdx));
        xc = X{n}(CptsNdx);
        zc = Z{n}(CptsNdx);
        j = find(zc<=maxZ);

        if ~isempty(j)
            vbb(1) = min(vbb(1), min(zc(j)));
            vbb(2) = max(vbb(2), max(zc(j)));
            vbb(3) = min(vbb(3), min(xc(j)));
            vbb(4) = max(vbb(4), max(xc(j)));
        end
    end
end
%keyboard

% Fit volume so that ground fits in an image of size [1024*780] and rescale
% all objects
Dx = -vbb(3);
Dz = -vbb(1);
Dy = 0;
S = (booksize-1)/max(vbb(2)-vbb(1), vbb(4)-vbb(3));
for n = 1:Nobjects
    X{n} = fix(S*(X{n}+Dx))+1;
    Y{n} = fix(S*(Y{n}+Dy))+1;
    Z{n} = fix(S*(Z{n}+Dz))+1;
end
vbb = fix(S*[vbb(1)+Dz vbb(2)+Dz vbb(3)+Dx vbb(4)+Dx])+1;

% Create ground image
ground = 255*ones([vbb(2) vbb(4) 3], 'uint8');
for n = 1:Nobjects
% $$$     if strcmp(annotation.object(n).mesh.type, 'ground')
    if strcmp(annotation.object(n).polygon.polyType, 'ground')
        disp(annotation.object(n).name)
        ground = imageWarp(TX{n}*size(textures{n},2), TY{n}*size(textures{n},1), textures{n}, X{n}, Z{n}, ground, boundary{n});
    end
end

% Create standing models
m = 0;
clear book
objectColors = uint8(255*hsv(Nobjects));
for n = 1:Nobjects
  n
% $$$     if strcmp(annotation.object(n).mesh.type, 'foldedcard');
    if strcmp(annotation.object(n).polygon.polyType, 'standing');
        CptsNdx = boundary{n}(find(annotation.object(n).mesh.CptsNdx));
        txc = TX{n}(CptsNdx)*size(textures{n},2);
        xc = X{n}(CptsNdx);
        zc = Z{n}(CptsNdx);

        if max(X{n})-min(X{n})>5 & median(Z{n})<(S*maxZ+Dz) % only render objects with at least 5 pixels width
            % sort contact points
            [txc,v] = sort(txc);
            xc = xc(v);
            zc = zc(v);
            
            if length(xc)==1
                % it is a fronto parallel object already
                Yp = -Y{n};
                Xp = X{n};
                w = max(5,(max(Xp)-min(Xp))/10);
                xc_p = [xc-w xc+w];
                xc = [xc-w xc+w];
                zc = [zc zc];
            else
                % Get frontal view for contact points
                D  = sqrt((xc(2:end)-xc(1:end-1)).^2 + (zc(2:end)-zc(1:end-1)).^2);
                xc_p = cumsum([1 D]);

                % Now warp the mesh
                %dx = abs(xc(2:end)-xc(1:end-1));
                %u = unique([1 find(dx>0) length(xc)]);
                %Xp = interp1(xc(u), xc_p(u), X{n}, 'linear', 'extrap');
                %Xp = interp1(xc, xc_p, X{n}, 'linear', 'extrap');
                Xp = interp1(txc, xc_p, TX{n}*size(textures{n},2), 'linear', 'extrap'); % This is probably wrong
                Yp = -Y{n};
            end
            
            
            xc_p = xc_p - min(Xp)+1;
            Yp = Yp-min(Yp)+1;
            Xp = Xp-min(Xp)+1;

            RENDER_OBJ = 1;%0;
            if RENDER_OBJ
              % warp the image
              object = 255*ones([ceil(max(Yp)) ceil(max(Xp)) 3], 'uint8');
              object = imageWarp(TX{n}*size(textures{n},2), TY{n}*size(textures{n},1), textures{n}, Xp, Yp, object, boundary{n});
            end
            
            % Make markings on the ground plane (interpolate between
            % consecutive contact points)
            for k = 1:length(xc)-1
                ground = drawline(zc(k),xc(k),zc(k+1),xc(k+1),ground,objectColors(n,:),3);
            end

            if RENDER_OBJ
              % Add end to glue together
              foot = getFoot(xc_p, size(object,2), booksize, objectColors(n,:));
              
              m = m+1;
              book{m} = cat(1,object,foot);
            end
            
            %figure
            %imshow(book{m})
            %pause
        end
    end
end
ground = ground(end:-1:1,:,:);
ground(:,1,:) = 0;
ground(1,:,:) = 0;
ground(:,end,:) = 0;
ground(end,:,:) = 0;
book{m+1} = ground;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compose book
N = length(book);
s = [];
for i = 1:N
    s(i,:) = size(book{i});
end


N = length(book);
nx = ceil(sqrt(N));
ny = ceil(N/nx);

areaBest = Inf;
for iter = 1:20
    BOOK = [];
    k = 0;
    ndx = randperm(N);
    for n = 1:ny
        k = k+1;
        if k<=N
            L = book{ndx(k)};
        end
        for m = 2:nx
            k = k+1;
            if k<=N
                L = addleft(L, 255*ones([5 5 3]));
                L = addleft(L, book{ndx(k)});
            end
        end
        if n == 1
            BOOK = L;
        else
            BOOK = addbottom(BOOK,255*ones([5 5 3]));
            BOOK = addbottom(BOOK,L);
        end
    end
    area = max(size(BOOK,1),size(BOOK,2)*1.2);
    if area < areaBest
        BOOKbest = BOOK;
        areaBest = area;
    end
end
BOOK = BOOKbest;

addHeader = 0;
if addHeader
  % Add header
  labelmeIcon = imread('LabelMeNEWtight198x55.jpg');
  HOMEIMAGES = 'http://labelme.csail.mit.edu/Images';
  thumb = LMsceneThumbnail(annotation, HOMEIMAGES, img, 'lines');
  ncb = size(BOOK,2);
  nci = size(labelmeIcon,2);
  nct = size(thumb,2);
  labelmeIcon = addleft(labelmeIcon, 255*ones([5 ncb-nci-nct-1 3]));
  labelmeIcon = addleft(labelmeIcon, thumb);
  labelmeIcon = addbottom(labelmeIcon, 255*ones([5 5 3]));
  labelmeIcon = addbottom(labelmeIcon, 128*ones([5 size(labelmeIcon,2) 3]));
  labelmeIcon = addbottom(labelmeIcon, 255*ones([5 5 3]));
  BOOK = addbottom(labelmeIcon, BOOK);
end

function img = drawline(x1,y1,x2,y2,img,color,thikness)
%
% create an line in the image

[nrows ncols cc] = size(img);

%cB = strel('disk',thikness);
cB = double((hamming(1+2*thikness)*hamming(1+2*thikness)')>.1);
m = ceil(max(abs(x2-x1), abs(y2-y1)));

x = round(linspace(x1,x2,m));
y = round(linspace(y1,y2,m));

% remove points outside the image
v = find(x>0 & y>0 & x<=nrows & y<=ncols);
x = x(v); y = y(v);

% draw line
mask = zeros([size(img,1) size(img,2)], 'single');
j = sub2ind(size(mask), x,y);
mask(j) = 1;

%mask = uint8(imdilate(mask, cB)==0);
mask = uint8(conv2(mask, cB, 'same')>0.1);


cmask = zeros([size(mask,1) size(mask,2) 3], 'uint8');
cmask(:,:,1) = (mask)*color(1);
cmask(:,:,2) = (mask)*color(2);
cmask(:,:,3) = (mask)*color(3);

mask = repmat(1-mask,[1 1 3]);
img = img.*mask+cmask;


function foot = getFoot(xc_p, ncols, booksize, linecolor)

s = 30*booksize/640;
foot = 255*ones([s+1 ncols 3], 'uint8');

Nc = length(xc_p);
for k = 1:Nc-1
    dx = (xc_p(k+1)-xc_p(k))/5;
    foot = drawline(1,xc_p(k),s,xc_p(k)+dx,foot,linecolor,1);
    foot = drawline(s,xc_p(k)+dx,s,xc_p(k+1)-dx,foot,linecolor,1);
    foot = drawline(1,xc_p(k+1),s,xc_p(k+1)-dx,foot,linecolor,1);
end



function img = addleft(img1,img2)

[nr1, nc1, c1] = size(img1);
[nr2, nc2, c2] = size(img2);

if nr1>nr2
    img2 = cat(1,img2, 255*ones([nr1-nr2 nc2 3]));
else
    img1 = cat(1,img1, 255*ones([nr2-nr1 nc1 3]));
end

img = cat(2,img1,img2);


function img = addbottom(img1,img2)

[nr1, nc1, c1] = size(img1);
[nr2, nc2, c2] = size(img2);

if nc1>nc2
    img2 = cat(2,img2, 255*ones([nr2 nc1-nc2 3]));
else
    img1 = cat(2,img1, 255*ones([nr1 nc2-nc1 3]));
end

img = cat(1,img1,img2);







