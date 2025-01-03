function [Q, R] = generateRotationMatrix(L)
%Generates Rotation Matrix out of the cholesky decomposition of the correlation matrix by rotating the
%bisector into e_1 via a basis-transformation matrix. 
N = bisector(L);
K = length(N);
E = triu(ones(K,K));
T = [N, E(:, 2:K)];
Q = gs(T);
Q = [Q(:, 2:K), Q(:,1)];
if Q(:, K) ~= N
    Q = -Q;
end
R = L * Q;
R_T = normc(inv(R));
assert(prod(abs(R_T(end, :) - R_T(end, 1)) < 0.00001), 'RoationError:ResultingMatrixNotParallel')
end


