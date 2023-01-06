function [U, J] = SRRQRtol(A, reltol)

QR_tol = reltol * max(vecnorm(A, 2));
[R, p, r] = QRpartial(A', 'tol', QR_tol);
R = R(1 : r, :);
J = p(1 : r);
P = [];
P(p) = 1:length(p);

[J, I] = sort(J, 'ascend');
invI = zeros(length(I),1);
invI(I) = 1 : length(I);
invI = [invI', (length(I)+1):length(P)];
P = invI(P);

E = linsolve(R(1:r, 1:r), R(1:r, (r+1):end), struct('UT', true))';
E = E(:, I);
U = [eye(size(E, 2)); E];
U = U(P, :);

end