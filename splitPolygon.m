function [Txn, Tyn, Xn, Yn, Zn] = splitPolygon(Tx,Ty,X,Y,Z)
%
 
%keyboard

Tx = [Tx; Tx(1)];
Ty = [Ty; Ty(1)];
X = [X; X(1)];
Y = [Y; Y(1)];
Z = [Z; Z(1)];



[h,xnn]=hist(Tx,unique(Tx));
j = find(h==2);
xc = xnn(j);


Xn{1} = X;
Yn{1} = Y;
Zn{1} = Z;
Txn{1} = Tx;
Tyn{1} = Ty;
nc = length(xc);

for m = 1:nc
    j1 = find(Txn{m}<=xc(m));
    j2 = find(Txn{m}>=xc(m));
    
    %keyboard
    
    Xc = Xn{m};
    Yc = Yn{m};
    Zc = Zn{m};
    Txc = Txn{m};
    Tyc = Tyn{m};
    
    Xn{m} = Xc(j1);
    Yn{m} = Yc(j1);
    Zn{m} = Zc(j1);
    Txn{m} = Txc(j1);
    Tyn{m} = Tyc(j1);
    
    Xn{m+1} = Xc(j2);
    Yn{m+1} = Yc(j2);
    Zn{m+1} = Zc(j2);
    Txn{m+1} = Txc(j2);
    Tyn{m+1} = Tyc(j2);
end


if 0
    figure
    subplot(221)
    plot3(X,Y,Z,'o-')
    subplot(222)
    for m = 1:nc
        plot3(Xn{m},Yn{m},Zn{m},'o-')
        hold on
    end
    subplot(223)
    plot(Tx,Ty,'o-')
    subplot(224)
    for m = 1:nc
        plot(Txn{m},Tyn{m},'o-')
        hold on
    end
end
