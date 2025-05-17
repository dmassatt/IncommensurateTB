% Produces a matrix operating on the dof space (1:N)^2 corresponding to a
% shift in x up, i.e. psi at x = n is moved to the x = n+1 slot.

function Tx = Generate_Tx(N)

Index = 1:N^2;
Y_shift = mod(Index + N-1, N^2)+1;

Tx = sparse(Index,Y_shift,ones(size(Index)),N^2,N^2);

end