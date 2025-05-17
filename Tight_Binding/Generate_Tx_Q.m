% Produces a matrix operating on the dof space (1:N)^2 corresponding to a
% shift in x up, i.e. psi at x = n is moved to the x = n+1 slot.

% This is a shift with momenta Q Bloch boundary conditions.

function Tx = Generate_Tx_Q(N, Q)

Index = 1:N^2;
Y_shift = mod(Index + N-1, N^2)+1;
S = [ones(1,N^2-N), ones(1,N)*exp(1i*Q)];


Tx = sparse(Index,Y_shift,S,N^2,N^2);

end