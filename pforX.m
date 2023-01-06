function [ mpd, Q ]  = pforX( hX, p )
% determine adaptive tensor structure
% INPUT:
%       hX = [hx1,hx2,...] length of each dimension of X
%       p: approximation level  (p>=0)
% OUTPUT:
%       mpd = [p1,...,pd] with p1+...+pd = p
%       Q = (1+p1)*...*(1+pd)   total #nodes

pX = p*hX/sum(hX);
mpd = floor(pX);

while sum(mpd) < p
    [~,ind] = max(pX-mpd);
    if length(ind) == 1
        mpd(ind) = mpd(ind) + 1;
    else
        for k = 1:min( length(ind), p-sum(mpd) )
            mpd(ind(k)) = mpd(ind(k)) + 1;
        end
    end
end
Q = prod(1+mpd);
