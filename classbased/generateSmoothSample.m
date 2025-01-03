function [samples, probs] = generateSmoothSample(L, B, U)
%input: L: KxK invertable Lower triangular matrix
%B: Kx1 vector
%U: 1xK vector of random UNIF[0,1] numbers
%output: samples: KX1 Vector
%probs: 1xK vector
%Generates a sample, and the corresponding likelihood from a polyhedral as
%defined by Lx>=B with L being a triangular matrix. 
K = length(L);
samples = zeros(K, 1);
probs = ones(1, K);
for k = 1:K
    B(k)  =  B(k) - L(k,1:k-1) * samples(1 : k - 1);
    if L(k,k) < 0
        lbound = - inf;
        ubound = B(k) / L(k,k);
    else
        lbound = B(k) / L(k,k);
        ubound = inf;
    end 
    samples(k) = truncNormSample(lbound, ubound, U(k));
    probs(k) = probNormSlice(lbound, ubound);
end
end

