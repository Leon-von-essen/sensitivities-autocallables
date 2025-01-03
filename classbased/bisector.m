function vec = bisector(L)
%Inputs:
%L: KXK lower triangular matrix with L_11 = 1
%returns: 
%vec: KX1 vector
%Generates a normalized bisector vector vec from an invertable Matrix L, that describes
%a convex Cone via LX >= 0. First we get the edges as the columns of the
%inverse of L
%(v_i = inv(L)*e_i), as this set K-1 equalities to 0, while leaving one
%open, corresponding to the intersection of K-1 K-1-dimensional hyperplanes.
%with V = inv(L) - v_1 and x = inv(V) * e_1, we ensure Lx >= 0 with V(1,:) = e_1, as
%L_11>0.
K = length(L);
edges = normc(inv(L))';
V = edges - edges(1,:);
res = zeros(1, K);
res(K) = 1;
V(1,:) = res;
V = [V(2:end, :); V(1, :)];
n = V \ res';
vec = n / norm(n);
end


