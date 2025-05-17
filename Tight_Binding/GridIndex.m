% returns the index of the meshgrid corresponding to data (x,y) on grid
% (1:N)^2.

function index = GridIndex(x,y,N)

index = N*(x-1) + y;

end