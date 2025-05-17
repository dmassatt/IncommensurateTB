% Produces a matrix operating on the dof space (1:N)^2 corresponding to a
% shift in y to the right, i.e. psi at y = n is moved to the y = n+1 slot.

% This is a shift with momenta Q Bloch boundary conditions.

function Ty = Generate_Ty_Q(N,Q)

Index = 1:N^2;
X_shift = mod(Index , N) - mod(Index-1, N) + Index;

S = ones(1,N^2);
S(N:N:(N^2)) = S(N:N:(N^2)) * exp(1i*Q);

Ty = sparse(Index,X_shift,S,N^2,N^2);

end