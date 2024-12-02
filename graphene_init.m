function loc_s = graphene_init(theta,f1,f2,r_cut)

L = [1 cos(pi/3);0 sin(pi/3)]*(1.43*sqrt(3));

basis = [0 1; 0 1/sqrt(3)]*(1.43*sqrt(3));

R = [cos(theta) -sin(theta);sin(theta) cos(theta)];

M = ceil(2*r_cut/(norm(L(:,1))*sin(pi/3)));
N = ceil(2*r_cut/(norm(L(:,2))*sin(pi/3)));

if mod(M,2)==0
    M = M+1;
end
if mod(N,2)==0
    N = N+1;
end
Size = [M,N];
s1 = sheet(L,Size,basis,r_cut,[0;0]);
s2 = sheet(R*L,Size,R*basis,r_cut,[0;0]);
intra_c = 4.33;
inter_c = 6 + norm(L*[1;1]);

intra = @(x,y,o,m) intra_graphene(x,y,o,m);
%inter = @(x,y,o,m,a,b) inter_graphene(x,y,o,m,a,b,theta);
%s1.view('b');
% hold on
% s2.view('g');
% scatter(r_cut*cos((1:1:100)*pi/50),r_cut*sin((1:1:100)*pi/50),'r')
% hold off

loc_s = loc_system(s1,s2,intra,intra,theta,intra_c,intra_c,inter_c,f1,f2);
end

