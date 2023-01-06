function [ U, V ] = geoLR( X, Y, pk, ff )
%@@@ Fast low-rank approximation of general kernel matrices
%===    "Data-Driven Linear Complexity Low-Rank Approximation 
%===    of General Kernel Matrices: A Geometric Approach"
%===        by D.Cai, E.Chow, Y.Xi 
%===         arXiv:2212.12674
%--------------------------------------------------------
% INPUT:
%       X: n points in |R^d -> n-by-d matrix
%       p: approximation level (not rank!)
%       ff: kernel function handle, say ff = @(x,y) pdist2(x,y).^2
% OUTPUT:
%       A ~~~ U*V

[ ind, ~ ] = AnchorNetApp( Y, pk, 'U', ff );
r_ANC = length(ind);
tol = 1e-13; % tunable
[U,jdex] = SRRQRtol(ff(X,Y(ind,:)), tol);
V = ff(X(jdex,:),Y);
r_out = length(jdex);
fprintf('app_level=%d, Out rank=%d, ANC rank=%d\n', pk,r_out,r_ANC)