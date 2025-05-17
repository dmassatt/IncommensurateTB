% returns Bilayer Graphene spatial position vector, R + \tau, for all dofs.

function [X,Y] = PositionVecs(L)

[Nx,Ny] = meshgrid(-L:L,-L:L);   % grid

A = [3/2  3/2; sqrt(3)/2 -sqrt(3)/2];

XA = A(1,1)*Nx + A(1,2)*Ny;
YA = A(2,1)*Nx + A(2,2)*Ny;
Uno = ones(size(XA));
XB = XA + Uno;
YB = YA;
Uno = Uno(:);

X = [XA(:);XB(:);XA(:)+Uno; XB(:)+Uno];
Y = [YA(:);YB(:);YA(:);YB(:)];

end