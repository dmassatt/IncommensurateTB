function P = Bilayer_Projector(Sheet1, Sheet2, f1,f2,shift)
% f = f([x;y])
if ~exist('shift','var')
    shift = [0;0];
end
Size1 = size(Sheet1.Atom_Positions,2);
Size2 = size(Sheet2.Atom_Positions,2);
Size = Size1 + Size2;

P = zeros(Size,1);
%Smon.Basis_Size

for m = 1:Sheet1.Sheet_Size(1,1)
    for n = 1:Sheet1.Sheet_Size(1,2)
        for k = 1:Sheet1.Basis_Size
            if f1(Sheet1.Atom_Positions(:,Sheet1.index(m,n,k))-Sheet1.Origin+shift)
                P(Sheet1.index(m,n,k),1) = 1;
            end
        end
    end
end
for m = 1:Sheet2.Sheet_Size(1,1)
    for n = 1:Sheet2.Sheet_Size(1,2)
        for k = 1:Sheet2.Basis_Size
            if f2(Sheet2.Atom_Positions(:,Sheet2.index(m,n,k))-Sheet2.Origin)
                P(Size1 + Sheet2.index(m,n,k),1) = 1;
            end
        end
    end
end
end