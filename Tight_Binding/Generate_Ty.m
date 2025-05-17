% Produces a matrix operating on the dof space (1:N)^2 corresponding to a
% shift in y to the right, i.e. psi at y = n is moved to the y = n+1 slot.

function Ty = Generate_Ty(N)

Index = 1:N^2;
X_shift = mod(Index , N) - mod(Index-1, N) + Index; 

Ty = sparse(Index,X_shift,ones(size(Index)),N^2,N^2);

end