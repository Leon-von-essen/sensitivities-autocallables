function [T, R, proj_Rs, proj_bs, Qs] = loadProjectedMatrices(L)
%%Input: 
%L: KxK cholesky decomposition of the correlation matrix
%Output:
%R: KxK matrix without zeros that has RX ~ N(0, correlationMatrix)
%volatilities: Kx1 vector with the stddevs of the covMatrix on the diag
%proj_Rs: K^2-KxK^2-K matrix of projected blockmatrices 
%proj_bs: K^2-KxK matrix of boundary condition projection blockmatrices.
%Qs: K^2-KxK^2-K matrix of rotation submatrices (unitary) 
%Loads projected matrices from a covariance matrix
[T, R] = generateRotationMatrix(L);
[proj_Rs, proj_bs, Qs] = altBarrierProjections(R);
end