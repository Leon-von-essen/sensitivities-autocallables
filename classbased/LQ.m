function [L, Q] = LQ(A)
if ~all(all(A == tril(A)))
[U, R] = qr(A');
L = R';
Q = U';
else 
    L = A;
    Q = eye(size(A));
end
end
