
D = -5:.01:5;

orbital1 = 1;
orbital2 = 2; % 1 = A site, 2 = B site
[X,Y] = meshgrid(D,D);

h = inter_graphene(X,Y,orbital1,orbital2,6*pi/180);

surf(X,Y,h,'EdgeColor','none') % plots hopping function as function of (x,y)

view(2)