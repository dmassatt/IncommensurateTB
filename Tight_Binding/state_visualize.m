% given a state existing on a (2*L+1)^2 lattice, visualize a centered
% (2*W+1)^2 grid.
function value = state_visualize(L,W, psi,background)
if ~exist('background')
    background = 1e-8;
end

N = 2*L+1;
KeepGrid = (L+1 - W):(L+1 + W);

v1 = psi(1:N^2);
v2 = psi((N^2+1):(2*N^2));
v3 = psi((2*N^2+1):(3*N^2));
v4 = psi((3*N^2+1):(4*N^2));
v1 = reshape(v1,N,N);
v2 = reshape(v2,N,N);
v3 = reshape(v3,N,N);
v4 = reshape(v4,N,N);

V1 = v1(KeepGrid,KeepGrid);
V2 = v2(KeepGrid,KeepGrid);
V3 = v3(KeepGrid,KeepGrid);
V4 = v4(KeepGrid,KeepGrid);

[Nx,Ny] = meshgrid(-W:W,-W:W);   % grid

A = [3/2  3/2; sqrt(3)/2 -sqrt(3)/2];

XA = A(1,1)*Nx + A(1,2)*Ny;
YA = A(2,1)*Nx + A(2,2)*Ny;

XB = XA - ones(size(XA));
YB = YA;

zmax = max([abs(V1(:));abs(V2(:));abs(V3(:));abs(V4(:))]);

PointMaxSize = 200;
ZA = abs(PointMaxSize*V1/zmax);
ZB = abs(PointMaxSize*V2/zmax);
ZA2 = abs(PointMaxSize*V3/zmax);
ZB2 = abs(PointMaxSize*V4/zmax);
background = 1e-8;
hold off
figure(1)
scatter(XA(:),YA(:),ZA(:)+background,'fill','b')
hold on
scatter(XB(:), YB(:),ZB(:)+background,'fill','b')
scatter(-1,0,300,'r')
axis equal
title('layer 1')
hold off

figure(2)
scatter(XA(:),YA(:),ZA2(:)+background,'fill','r')
hold on
scatter(XB(:), YB(:),ZB2(:)+background,'fill','r')
scatter(0,0,300,'b')
axis equal
hold off
title('layer 2')

%figure(4)
%semilogy(-W:W,ZA(:,W+1))

end