%========================================================
%===  "Data-Driven Linear Complexity Low-Rank Approximation 
%===  of General Kernel Matrices: A Geometric Approach"
%===      by D.Cai, E.Chow, Y.Xi 
%===       arXiv:2212.12674
%========================================================

%-----  All points are in |R^d
d = 3; M = 2e3; N = 3e3;
X = randn(M,d); Y = randn(N,d); % row size = #pts

%----- assign a suitable bandwidth parameter
L = 0.3*max( pdist2(X(1:min(M,30),:),Y,'euclidean'), [], 'All' );
ff = @(x,y) exp(-pdist2(x,y,'euclidean').^2/L^2);
% Make sure the kernel matrix is numerically low-rank
% Otherwise it's meaningless to use low-rank approximation

%----- Test Data-driven Geomtric Low-Rank Compression
A = ff(X,Y); % for checking approximation error only

for pk = 1:5:30
    % pk controls the approximation level (shape of anchor net)
    % NOTE: pk is not equal to rank or #pts selected!

    
%>>>>>>>>>    geometric Low-Rank compression    <<<<<<<<<
    [ U, V ] = geoLR( X, Y, pk, ff );
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    tmp = (A-U*V).^2;
    Err_max = max( sqrt(tmp(:)) );
    Err_Fro = sqrt(sum(tmp(:)));
    fprintf('||abs_err||_max = %5.2e, ||abs_err||_F = %5.2e\n',Err_max,Err_Fro)
end

%----- FAQ: No accuracy when rank is large?
% Check singular values of A=ff(X,Y) to see
% whether they have rapid enough decay.
% If not, low-rank approximation won't be accurate.

