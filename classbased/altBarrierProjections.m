function [Ls, proj_bs, Qs] = altBarrierProjections(R)
%Input: 
%R: KxK nonzero matrix
%Output:
%Ls: K ^ 2xK ^ 2 lower triangular
%proj_bs: K ^ 2xK ^ 2
%Qs:K ^ 2xK ^ 2 Orthonormal
%Generates the projections that take a H-representation of a polyhedral Rx
%>= B and projects one boundary hyperplane
%onto the K-1 first coordinates. The boundary is identified by setting 
%the row-th equation to =. Then reenters the first boudnary condition
assert(all(R(:,1)), 'projectedConditions:FirstColumnNonZero', ...
       'First Column of matrix contains 0');
K = length(R);
proj_Rs = zeros(K ^ 2, K ^ 2);
proj_bs = zeros(K ^ 2, K ^ 2);
for row = 1:K
    if row == 1
    swap_m = rot90(eye(K, K));
    else 
        swap_m = eye(K, K);
    end
    proj_R = R - 1 / R(row, K) * R(:, K) * R(row, :); 
    proj_R(row, :)  = - R(row, :);
    proj_R = swap_m * proj_R;
    proj_b = eye(K) - [zeros(K, row - 1), R(:, K), ...
                       zeros(K, K - row)] / R(row, K);
    proj_b(row, :) = [zeros(1, row - 1), -1, zeros(1, K - row)];
    proj_b = swap_m * proj_b; 
    row_ind = (row - 1) * K + 1 : row * K;
    proj_bs(row_ind, row_ind) = proj_b;
    proj_Rs(row_ind, row_ind) = proj_R;
end
[Ls, Qs] = LQ(proj_Rs);
end
