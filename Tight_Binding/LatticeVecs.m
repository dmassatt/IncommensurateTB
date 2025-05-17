function [X,Y] = LatticeVecs(L)

[Nx,Ny] = meshgrid(-L:L,-L:L);   % grid

A = [3/2  3/2; sqrt(3)/2 -sqrt(3)/2];

X = A(1,1)*Nx + A(1,2)*Ny;
Y = A(2,1)*Nx + A(2,2)*Ny;

end