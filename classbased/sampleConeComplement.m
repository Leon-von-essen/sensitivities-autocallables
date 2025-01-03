function [samples, probs] = sampleConeComplement(proj_Rs, proj_bs, Qs,...
                                                 B, U)
%Input: 
%R: KxK Matrix without zeros
%B: Kx1-Vector 
%Output:
%samples: KxK matrix of truncated samples
%probs: KxK matrix of corresponding probabilities
%U:1xK vector of UNIF[0,1] variables (iid). 
%Samples the complement of a convex cone in a smooth manner. The cone is
%given in its H-representation R >= B. 
%probs are sampled in inverse order and are reordered in the end for
%consistency
K = length(B);
samples = zeros(K, K);
probs = zeros(K, K);
for row = 1:K
    r_ind = (row - 1) * K + 1: row  * K;
    R_proj = proj_Rs(r_ind , r_ind);
    Q = Qs(r_ind , r_ind);
    B_proj = proj_bs(r_ind, r_ind) * B';
    
    [samples_pre, ...
     probs(row, :)]   = generateSmoothSample(R_proj, B_proj, U(:, row));
    assert(all(R_proj * samples_pre > B_proj), 'sampling_error: sample outside of barrier')
    samples(row, :) = samples_pre' * Q;
end
end