function [Psup,objectClasses] = supportGraph(DB, objectClasses, Ppart)
%
% Inputs:
%    DB - LabelMe DB structure.
%
% Outputs:
%    Psup(i,j) - Probability of object class Oi supports Oi.
%    objectClasses - Cell array of object classes corresponding to the
%                    rows/cols of Psup.

% Parameters:
thresh_par = 0.5; % parts threshold (to collect local evidences)
thresh_sup = 25; % support threshold (to collect local evidences)
alpha = 5; % pseudocounts

% Only use images that are mostly labeled:
%DB = DB(LMlabeledarea(DB)>=thresh_labeled_area);

% Extract object classes:
if nargin == 1
    objectClasses = LMobjectnames(DB,'name');
end
objectClasses = lower(objectClasses);
Nobjects = length(objectClasses);
Nimages = length(DB);

% Compute part probabilities
PC_s = zeros([Nobjects, Nobjects]); % P(overlap | no part)
Counts = zeros([Nobjects, Nobjects]); % P(overlap | no part)
for i = 1:Nimages
    display(sprintf('%d out of %d',i,Nimages));

    %keyboard
    
    if isfield(DB(i).annotation,'object')
        % Get object indices for this image
        Nobj = length(DB(i).annotation.object);
        clear ndx X Y x y
        for j = 1:Nobj
            ndx(j) = strmatch(strtrim(lower(DB(i).annotation.object(j).name)), objectClasses, 'exact');
            [X{j},Y{j}] = getLMpolygon(DB(i).annotation.object(j).polygon);
            [x(j),y(j)] = GetBottomPoint(X{j},Y{j});
        end

        % Get Pso for this image
        for m = 1:Nobj
            for cl = unique(ndx) % loop on classes in this image (each class counts only once as background)
                if Ppart(cl,ndx(m))<.9
                    n = find(ndx==cl);

                    clear C
                    for k = 1:length(n)
                        C(k) = inpolygon(x(m),y(m), X{n(k)},Y{n(k)});
                    end
                    C = sum(C)>0;
                    
                    n = n(1);
                    PC_s(ndx(n),ndx(m)) = PC_s(ndx(n),ndx(m))+C;
                    Counts(ndx(n),ndx(m)) = Counts(ndx(n),ndx(m))+1;
                end
            end
        end
    end
end

Psup = PC_s./(Counts+eps);
Psup = Psup .* (Counts>5);
Psup = Psup - diag(diag(Psup));

% Just for visualization:
[i,j] = find(Psup>.5 & Counts>10); %only show results when there are enough samples
v = Psup(sub2ind(size(Psup),i,j));
[v,n] = sort(v, 'descend');
for k=n'
    disp(sprintf('%s supports %s (with probability = %1.2f, estimated from %d samples)', objectClasses{i(k)}, objectClasses{j(k)}, Psup(i(k),j(k)), Counts(i(k),j(k))))
end
