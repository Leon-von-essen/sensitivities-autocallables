function Q = gs(T)
T = normc(T);
K = length(T);
Q = zeros(K, K);
Q(:, 1) = T(:, 1);
for k = 2:K
    Q(:, k) = T(:, k);
    for j = 1:k-1
        Q(:, k) = Q(:, k) - (Q(:, j)' * T(:, k)) * Q(:, j);
    end
    Q = normc(Q);
end

end