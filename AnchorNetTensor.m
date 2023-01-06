function [ Xref ] = AnchorNetTensor( X, p, choice, ver )
% Given point set X and approximation level p >= 0,
% generate Q anchor points Xref
% INPUT:
%       X: n-by-d matrix
%       p: approximation level 0,1,2,...
%       choice: type of anchors
% OUTPUT:
%       Xref: anchor net pts (Q-by-d)

[n,d] = size(X);

Xmin = min(X); Xmax = max(X);
hX = Xmax-Xmin; %  hX(k) = 0 is unlikely to be true

[ mpd, Q ]  = pforX( hX, p );
pmod = p;
if ver == 1
    Qm = max(round(n),5e3);
else
    Qm = max(round(n/2),1e4);
end
while Q > Qm & n > 1e3 & d > 3
    pmod = pmod - 1;
    [ mpd, Q ]  = pforX( hX, pmod );
end




% Xrefd{i}=anchor pts in the i^th dimension
Xrefd = cell(1,d);
if min(mpd) == 0
    Xrefd(mpd==0) = num2cell( (Xmin(mpd==0)+Xmax(mpd==0))/2 );
end
ind = 1:d;
ind = ind(mpd>0);
%@@@@ choice=1,2,3...: anchors pts types
% larger choice -> more near the boundary
switch choice
    case 'U'
        for k = ind
            hk = hX(k)/mpd(k);
            Xrefd{k} = Xmin(k)+hk*(0:mpd(k))';
        end
    case 'C'
        for k = ind
            Xrefd{k} = (Xmin(k)+Xmax(k))/2+(Xmax(k)-Xmin(k))/2*cos( (2*(0:mpd(k))'+1)*pi/(2*mpd(k)+2) );
        end
    case 1
        for k = ind
            hk = hX(k)/(mpd(k)+1);
            Xrefd{k} = Xmin(k)+hk/2+hk*(0:mpd(k))'; % column vector !
        end
    case 2
        for k = ind
            hk = hX(k)/(mpd(k)+2/3);
            Xrefd{k} = Xmin(k)+1/3*hk+hk*(0:mpd(k))';
        end
    case 3
        for k = ind
            hk = hX(k)/(mpd(k)+1/2);
            Xrefd{k} = Xmin(k)+1/4*hk+hk*(0:mpd(k))';
        end
    case 4
        for k = ind
            hk = hX(k)/(mpd(k)+2/5);
            Xrefd{k} = Xmin(k)+1/5*hk+hk*(0:mpd(k))';
        end
    case 5
        for k = ind
            hk = 0.85*hX(k)/mpd(k);
            Xrefd{k} = Xmin(k)+0.1*hk+hk*(0:mpd(k))';
        end
end

% ready to build
Xref = zeros(Q,d);
if min(mpd) == 0
    tmp = cell2mat(Xrefd(mpd==0));
    Xref(:,mpd==0) = repmat(tmp,Q,1);
end
ind = 1:d;
ind = ind(mpd>0);
tk = Q;
for k = ind
    tk = tk/(1+mpd(k));
    Xref(:,k) = kron( ones(Q/(tk*(mpd(k)+1)),1), kron(Xrefd{k},ones(tk,1))  );
end
