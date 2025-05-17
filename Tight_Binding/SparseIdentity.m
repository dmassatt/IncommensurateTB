% generates a sparse matrix that is the identity on the N^2 degrees of
% freedom

function I = SparseIdentity(N)

J = 1:N;
I = J;
V = ones(1,N);

I = sparse(I,J,V);
end