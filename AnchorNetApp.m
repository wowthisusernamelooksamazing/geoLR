function [ xox, fac ] = AnchorNetApp( X, p, choice, fc )
%@@@ Fast low-rank approximation of symmetric kernel matrix
%===    "FAST DETERMINISTIC APPROXIMATION OF
%===    SYMMETRIC INDEFINITE KERNEL MATRICES WITH
%===    HIGH DIMENSIONAL DATASETS"
%===        by D.Cai, J.Nagy, Y.Xi
%--------------------------------------------------------
% INPUT:
%       X: n points in |R^d
%       p: approximation level 
%       choice: type of anchor points
%       fc: kernel function handle
% OUTPUT:
%       xox: indices of landmark points
%       fac{1}*fac{2} \approx A

[n,d] = size(X);
[Xmin, idmin] = min(X); [Xmax, idmax] = max(X);
hX = Xmax-Xmin;

if p < 3
    [~,i0] = min(Xmin);
    idmin = idmin(i0);
    [~,i0] = max(Xmax);
    idmax = idmax(i0);
end

% determine anchor net shape parameters
[ ~, Q ]  = pforX( hX, p );
p1 = p;
if Q > 3e3
    p = round(p/2);
    if d > 5
        p = min([p,10]);
    end
end
q = round(p/2+2);

% first round
[ Xref ] = AnchorNetTensor( X, p1, choice, 1 );
[~,tmp] = pdist2(X,Xref,'euclidean','Smallest',1);
xox0 = unique([reshape(tmp,1,prod(size(tmp))),idmin,idmax]);


if d > 5
    mmm = 10*p; % mmm ~ p---40*p
    [~,yo] = max(hX); [~,i1] = max(X(:,yo)); [~,i2] = min(X(:,yo));
    AN_1 = X(i1,:); AN_2 = X(i2,:); AN_3 = mean(X);
    m1 = floor(mmm/3); m2 = m1; m3 = mmm-m1-m2;
    dd1 = pdist2(X,AN_1); dd2 = pdist2(X,AN_2); dd3 = pdist2(X,AN_3);
    [~,iii1] = sort(dd1); [~,iii2] = sort(dd2); [~,iii3] = sort(dd3);
    id1 = iii1(1:round(n/m1):end); id2 = iii2(1:round(n/m2):end); id3 = iii3(1:round(n/m3):end);
    ss = unique([id1;id2;id3]);
    Xref = X(ss,:);
end

[~,ind] = pdist2(Xref,X,'euclidean','Smallest',1);

% second round
xox1 = [];
AN = cell(1,max(ind));
for i = unique(ind)
    gb = 1:n;
    indXi = gb(ind==i);
    Xi = X(indXi,:);
    ni = length(indXi);
    hXi = zeros(1,d);
    for iter = 1:d
        hXi(iter) = max(Xi(:,iter))-min(Xi(:,iter));
    end
    qi = min([q,ceil(q*max(hXi)/min(hX)),ceil(d*ni^(1/d)-d)]);
    Xref = AnchorNetTensor( Xi, max(qi-1,0), choice, 2 );
    [~,tmp] = pdist2(Xi,Xref,'euclidean','Smallest',1);
    xox_i = unique(reshape(tmp,1,prod(size(tmp))));
    xox1 = [xox1 indXi(xox_i)];
end

% pass
xox = unique([xox0,xox1]);
fac = cell(1,2);
X1 = X(xox,:);
fac{1} = fc(X,X1);
fac{2} = fc(X1,X1)\fac{1}.';
