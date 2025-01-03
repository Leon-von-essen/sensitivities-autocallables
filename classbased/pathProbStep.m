function [P, S] = pathProbStep(R, proj_Rs, proj_bs, Qs)
K = length(R);
syms s [K, 1]
syms r t rkk t_0
syms sigma [K, 1] 
syms u [K, 1] 
syms B b [K, 1]
syms X H [K, 1] 
syms p [K, 1] 
syms F [K, K]
syms P [K, K]
for k= 1:K
    loc = (k-1) * K + 1:K*k;
    pR = proj_Rs(loc, loc);
    T = tril(pR, -1);
    pB = proj_bs(loc, loc);
    Q = Qs(loc, loc);
    f = normcdf(1 ./ diag(pR) .* (pB * B - T * X));
    g = (log(b ./ s) - (r - (sigma .^ 2) / 2) * t) / sqrt(t) ./ sigma;
    h = UNSAFE_norminv((diag(pR) < 1) .* p .* u + ...
                       (diag(pR) < 1) .* p + (1 - p) .* u);
    
    F(:, 1) = subs(f, [B, X], [g, zeros(2,1)]);
    H = subs(h, p, F(:, 1));
    for j = 2:K
        F(:, j) = subs(f, [B, X], [g, H]);
        H = subs(h, p, F(:, j));
    end
    P(:, k) = diag(F);
    S(:, k) = s .* exp((r - (sigma .^ 2) / 2) * t + ...
                 sqrt(t) * sigma .* R * Q * H);
end
S = simplifyFraction(simplify(S));
P = simplifyFraction(simplify(P));
end