function [H,v] = GenerateH(theta,r_cut)
%(s1,s2,intra1,intra2,theta,intra_c1,intra_c2,inter_c,f1,f2)

f = @(x) norm(x) < r_cut;


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
inter_cut = 6 + norm(L*[1;1]);

intra = @(x,y,o,m) intra_graphene(x,y,o,m);

intra1_buffer = max([norm(s1.Lattice(:,1)),norm(s1.Lattice(:,2)),norm(s1.Lattice(:,1)+s1.Lattice(:,2))]);
intra2_buffer = max([norm(s2.Lattice(:,1)),norm(s2.Lattice(:,2)),norm(s2.Lattice(:,1)+s2.Lattice(:,2))]);
inter_buffer  = max([intra1_buffer,intra2_buffer]);

s1_size = size(s1.Atom_Positions,2);

s2_size = size(s2.Atom_Positions,2);

count = 0;

P = Bilayer_Projector(s1,s2,f,f);
mat_size = sum(P(:));

% Now it builds the first matrix
%% Counting
for i = 1:s1_size
    if P(i,1) == 1
        count = count + size(s1.neighbors(s1.Atom_Positions(:,i)-s1.Origin, intra_c+ intra1_buffer),1);
        count = count + size(s2.neighbors(s1.Atom_Positions(:,i)-s1.Origin, inter_cut+ inter_buffer),1);
    end
end
for i = 1:s2_size
    if P(i+s1_size,1) == 1
        count = count + size(s1.neighbors(s2.Atom_Positions(:,i)-s2.Origin,inter_cut+inter_buffer),1);
        count = count + size(s2.neighbors(s2.Atom_Positions(:,i)-s2.Origin,intra_c+intra2_buffer),1);
    end
end
%% Indexing
index = zeros(mat_size,2);
index_inv = zeros(max([s1_size,s2_size]),2);
index0 = 1;
for j = 1:(s1_size+s2_size)
    if  P(j,1) == 1
        f = j;          % sort out layers
        b = 1;
        if j > s1_size
            f = j - s1_size;
            b = 2;
        end
        index(index0,1) = f;
        index(index0,2) = b;
        index_inv(f,b) = index0;
        index0 = index0+1;
    end
end
I = zeros(count,1);
J = zeros(count,1);
S = J;
count2 = 1;

%% Calculate Matrix
% Calculate Cross-Terms

for i = 1:s2_size
    
    if  P(i+s1_size,1) == 1
        neighbors = s1.neighbors(s2.Atom_Positions(:,i)-s2.Origin, inter_cut+ inter_buffer);
        s = size(neighbors,1);
        X = count2:(count2+s-1);
        OneX = ones(size(X));
        I(X) =  index_inv(i,2)*OneX;
        J(X) =  index_inv(neighbors,1);
        
        r = (s2.Atom_Positions(:,i)-s2.Origin+s1.Origin)*OneX-s1.Atom_Positions(:,neighbors);

        S(X) =  interlayer_graphene(r(1,:),r(2,:),s2.index_inv(i,3)*OneX,s1.index_inv(neighbors,3)',2*OneX,OneX,theta);
        
        count2 = count2+s;
    end
end

for i = 1:s1_size
   
    if  P(i,1) == 1
        neighbors = s2.neighbors(s1.Atom_Positions(:,i)-s1.Origin, inter_cut+ inter_buffer);

        s = size(neighbors,1);
        X = count2:(count2+s-1);
        OneX = ones(size(X));
        I(X) =  index_inv(i,1)*OneX;
        J(X) =  index_inv(neighbors,2);
        r = (s1.Atom_Positions(:,i)-s1.Origin+s2.Origin)*OneX-s2.Atom_Positions(:,neighbors);
        
        S(X) =  interlayer_graphene(r(1,:),r(2,:),s1.index_inv(i,3)*OneX,s2.index_inv(neighbors,3)',OneX,2*OneX,theta);
        
        count2 = count2+s;

    end
end

inter_end = count2-1;
% Calculate Self Interaction

%S = zeros(size(I,1)- inter_end,1);

for i = 1:s1_size
    if  P(i,1) == 1
        neighbors = s1.neighbors(s1.Atom_Positions(:,i)-s1.Origin, intra_c+ intra1_buffer);
        s = size(neighbors,1);
        X = count2:(count2+s-1);
        OneX = ones(size(X));
        I(X) = index_inv(i,1)*OneX;
        J(X) = index_inv(neighbors,1);
        r = (s1.Atom_Positions(:,i)-s1.Origin+s1.Origin)*OneX-s1.Atom_Positions(:,neighbors);

        S(X) =  intra(r(1,:),r(2,:),s1.index_inv(i,3)*OneX,s1.index_inv(neighbors,3)');
        count2 = count2+s;

    end

end

for i = 1:s2_size
    if  P(i+s1_size,1) == 1
        neighbors = s2.neighbors(s2.Atom_Positions(:,i)-s2.Origin, intra_c+ intra2_buffer);
        s = size(neighbors,1);
        X = count2:(count2+s-1);
        OneX = ones(size(X));
        I(X) =  index_inv(i,2);
        J(X) =  index_inv(neighbors,2);
        r = (s2.Atom_Positions(:,i)-s2.Origin+s2.Origin)*OneX-s2.Atom_Positions(:,neighbors);

        S(X) =  intra(r(1,:),r(2,:),s2.index_inv(i,3)*OneX,s2.index_inv(neighbors,3)');
        
        count2 = count2+s;
    end
end

H = sparse(I,J,S,mat_size,mat_size);

v = zeros(2,2);
for s_n = 1:2
    for k = 1:2
        v(s_n,k) = index_inv(s1.Origin_Index(k,1),s_n);
    end
end

end