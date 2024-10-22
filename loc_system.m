%
% because of speed issues, I've built in the graphene interlayer coupling
% automatically. Function pointers are slow.

classdef loc_system
    properties
        sheet1
        sheet2
        intra_func1     % intra_func(x,y,o,m), (x,y) horizontal displacement, orbital o and m
        intra_func2
        %inter_func      % inter_func(x,y,o,m,a,b), (x,y) horizontal displacement vector,(o,m) are basis
        % index incase the atoms have different interaction
        % functions.
        % (a,b) are sheet numbers for (o,m) respectively
        % (x,y) positions for atom1-atom2
        
        intra1_cut      % cut-off radii for interactions
        intra2_cut
        intra1_buffer
        intra2_buffer
        inter_buffer
        inter_cut
        %R               % cut radius for region
        index           % these contain mappings between sheet indecies
        index_inv       % and matrix indecies
        
        % index(i,k) (i,1) = sheet index, (i,2) = sheet
        % index_inv(i,k) = matrix index (i = sheet index)
        % (k = sheet number)
        %H               % for zero-shift
        I               % sparsity pattern (I,J)
        J
        R1              % R1,R2 are x,y coordinates of R-R' distance between sites
        R2
        O1              % O1,O2 are orbitals that are interacting
        O2
        S1              % S1,S2 are sheets interacting
        S2
        %Sign            % Sign = 1 for sheet 1, Sign = -1 for sheet 2
        f1              % f1,f2 are total cut-off functions for sheet 1 & 2
        f2
        S               % INTRA layer entries only        
                        % sparsity value at [0,0]-shift
        P               % a vector of size of the parallelograms of atoms,
                        % defined such that P(k) = 1 if kth index included.
        inter_end
        
        theta
        
        mat_size
    end
    methods
        function obj = loc_system(s1,s2,intra1,intra2,theta,intra_c1,intra_c2,inter_c,f1,f2)

            obj.theta = theta;
            obj.sheet1      = s1;
            obj.sheet2      = s2;
            obj.f1          = f1;
            obj.f2          = f2;
            obj.intra_func1 = intra1;
            obj.intra_func2 = intra2;
            %obj.inter_func  = inter;
            obj.intra1_cut  = intra_c1;
            obj.intra2_cut  = intra_c2;
            obj.inter_cut   = inter_c;

            obj.intra1_buffer = max([norm(s1.Lattice(:,1)),norm(s1.Lattice(:,2)),norm(s1.Lattice(:,1)+s1.Lattice(:,2))]);
            obj.intra2_buffer = max([norm(s2.Lattice(:,1)),norm(s2.Lattice(:,2)),norm(s2.Lattice(:,1)+s2.Lattice(:,2))]);
            obj.inter_buffer  = max([obj.intra1_buffer,obj.intra2_buffer]);

            
            s1_size = size(s1.Atom_Positions,2);
            
            s2_size = size(s2.Atom_Positions,2);
            
            count = 0;
            
            obj.P = Bilayer_Projector(s1,s2,f1,f2);
            mat_size = sum(obj.P(:));
            obj.mat_size = mat_size;
            
            % Now it builds the first matrix
            %% Counting
            for i = 1:s1_size
                if obj.P(i,1) == 1
                    count = count + size(s1.neighbors(s1.Atom_Positions(:,i)-s1.Origin,obj.intra1_cut+obj.intra1_buffer),1);
                    count = count + size(s2.neighbors(s1.Atom_Positions(:,i)-s1.Origin,obj.inter_cut+obj.inter_buffer),1);
                end
            end
            for i = 1:s2_size
                if obj.P(i+s1_size,1) == 1
                    count = count + size(s1.neighbors(s2.Atom_Positions(:,i)-s2.Origin,obj.inter_cut+obj.inter_buffer),1);
                    count = count + size(s2.neighbors(s2.Atom_Positions(:,i)-s2.Origin,obj.intra2_cut+obj.intra2_buffer),1);
                end
            end
            %% Indexing
            obj.index = zeros(mat_size,2);
            obj.index_inv = zeros(max([s1_size,s2_size]),2);
            index0 = 1;
            for j = 1:(s1_size+s2_size)
                if obj.P(j,1) == 1
                    f = j;          % sort out layers
                    b = 1;
                    if j > s1_size
                        f = j - s1_size;
                        b = 2;
                    end
                    obj.index(index0,1) = f;
                    obj.index(index0,2) = b;
                    obj.index_inv(f,b) = index0;
                    index0 = index0+1;
                end
            end
            obj.I = zeros(count,1);%'int8');
            obj.J = zeros(count,1);%'int8');
            obj.S1 = zeros(count,1);%,'int8');
            obj.S2 = zeros(count,1);%,'int8');
            obj.R1 = zeros(count,1);
            obj.R2 = zeros(count,1);
            obj.O1 = zeros(count,1);%,'int8');
            obj.O2 = zeros(count,1);%,'int8');
            %obj.Sign = zeros(count,1);
            %obj.S = zeros(count,1);
            count2 = 1;
            
            %% Calculate Matrix
            % Then Calculate Cross-Terms
            
            for i = 1:s2_size
                
                if obj.P(i+s1_size,1) == 1
                    neighbors = s1.neighbors(s2.Atom_Positions(:,i)-s2.Origin,obj.inter_cut+obj.inter_buffer);
                    s = size(neighbors,1);
                    X = count2:(count2+s-1);
                    OneX = ones(size(X));
                    %for j = 1:size(neighbors,1)
                        obj.I(X) = obj.index_inv(i,2)*OneX;
                        obj.J(X) = obj.index_inv(neighbors,1);
                        
                        r = (s2.Atom_Positions(:,i)-s2.Origin+s1.Origin)*OneX-s1.Atom_Positions(:,neighbors);
                        obj.R1(X) = r(1,:);
                        obj.R2(X) = r(2,:);
                        obj.O1(X) = s2.index_inv(i,3);
                        obj.O2(X) = s1.index_inv(neighbors,3);
                        obj.S1(X) = 2;
                        obj.S2(X) = 1;
                        %obj.Sign(X) = -1;
                        
                        %obj.S(X) = obj.inter_func(r(1,:),r(2,:),s2.index_inv(i,3)*OneX,s1.index_inv(neighbors,3)',2*OneX,OneX);
                        
                        count2 = count2+s;
                    %end
                end
            end

            for i = 1:s1_size
               
                if obj.P(i,1) == 1
                    neighbors = s2.neighbors(s1.Atom_Positions(:,i)-s1.Origin,obj.inter_cut+obj.inter_buffer);

                    s = size(neighbors,1);
                    X = count2:(count2+s-1);
                    OneX = ones(size(X));
                    %for j = 1:size(neighbors,1)
                        obj.I(X) = obj.index_inv(i,1)*OneX;
                        obj.J(X) = obj.index_inv(neighbors,2);
                        r = (s1.Atom_Positions(:,i)-s1.Origin+s2.Origin)*OneX-s2.Atom_Positions(:,neighbors);
                        obj.R1(X) = r(1,:);
                        obj.R2(X) = r(2,:);
                        obj.O1(X) = s1.index_inv(i,3)*OneX;
                        obj.O2(X) = s2.index_inv(neighbors,3);
                        obj.S1(X) = OneX;
                        obj.S2(X) = 2*OneX;
                        %obj.Sign(X) = OneX;
                        
                        %obj.S(X) = obj.inter_func(r(1,:),r(2,:),s1.index_inv(i,3)*OneX,s2.index_inv(neighbors,3)',OneX,2*OneX);
                        
                        count2 = count2+s;
                    %end
                end
            end
 
            obj.inter_end = count2-1;
            % Calculate Self Interaction

            obj.S = zeros(size(obj.I,1)-obj.inter_end,1);
            
            for i = 1:s1_size
                if obj.P(i,1) == 1
                    neighbors = s1.neighbors(s1.Atom_Positions(:,i)-s1.Origin,obj.intra1_cut+obj.intra1_buffer);
                    s = size(neighbors,1);
                    X = count2:(count2+s-1);
                    OneX = ones(size(X));
                    %for j = 1:size(neighbors,1)
                        %if i ~= neighbors(j,1)
                            obj.I(X) = obj.index_inv(i,1)*OneX;
                            obj.J(X) = obj.index_inv(neighbors,1);
                            r = (s1.Atom_Positions(:,i)-s1.Origin+s1.Origin)*OneX-s1.Atom_Positions(:,neighbors);
                            obj.R1(X) = r(1,:);
                            obj.R2(X) = r(2,:);
                            obj.O1(X) = s1.index_inv(i,3)*OneX;
                            obj.O2(X) = s1.index_inv(neighbors,3);
                            obj.S(X-obj.inter_end) = obj.intra_func1(r(1,:),r(2,:),s1.index_inv(i,3)*OneX,s1.index_inv(neighbors,3)');
                            count2 = count2+s;
                        %end
                    %end
                end

            end

            for i = 1:s2_size
                if obj.P(i+s1_size,1) == 1
                    neighbors = s2.neighbors(s2.Atom_Positions(:,i)-s2.Origin,obj.intra2_cut+obj.intra2_buffer);
                    s = size(neighbors,1);
                    X = count2:(count2+s-1);
                    OneX = ones(size(X));
                    %for j = 1:size(neighbors,1)
                        %if i ~= neighbors(j,1)
                            obj.I(X) = obj.index_inv(i,2);
                            obj.J(X) = obj.index_inv(neighbors,2);
                            r = (s2.Atom_Positions(:,i)-s2.Origin+s2.Origin)*OneX-s2.Atom_Positions(:,neighbors);
                            obj.R1(X) = r(1,:);
                            obj.R2(X) = r(2,:);
                            obj.O1(X) = s2.index_inv(i,3)*OneX;
                            obj.O2(X) = s2.index_inv(neighbors,3);
                            obj.S(X-obj.inter_end) = obj.intra_func2(r(1,:),r(2,:),s2.index_inv(i,3)*OneX,s2.index_inv(neighbors,3)');
                            
                            count2 = count2+s;
                        %end
                    %end
                end
            end

            %[max(obj.I(:)), min(obj.I(:)) max(obj.J(:)), min(obj.J(:))]
            %obj.H = sparse(obj.I,obj.J,obj.S,mat_size,mat_size);

        end

        %%
        function X = PositionOperator(obj, j) % j = 1 or 2
                                                % generate position
                                                % operator for x,y
                                                % coordinate choice defined
                                                % by j.
            S = zeros(1,obj.mat_size);
            I = 1:obj.mat_size;
            for k = 1:obj.mat_size
                index = obj.index(k,:); % some position
                if index(2) == 2
                    v = obj.sheet2.Atom_Positions(:,index(1))-obj.sheet2.Origin;
                    S(k) = v(j);
                elseif index(2) == 1
                    v = obj.sheet1.Atom_Positions(:,index(1))-obj.sheet1.Origin;
                    S(k) = v(j);
                end
            end
            X = sparse(I,I,S,obj.mat_size,obj.mat_size);
        end
        
        function H_mat = MatrixShift(obj,shift)
            
            Shifts = shift*double(obj.S2-obj.S1)';
            r = [obj.R1';obj.R2']+Shifts;
            Sn = zeros(size(obj.I));
            Sn(1:obj.inter_end) =   obj.inter_func(r(1,1:obj.inter_end)',r(2,1:obj.inter_end)',obj.O1(1:obj.inter_end),...
                                    obj.O2(1:obj.inter_end),obj.S1(1:obj.inter_end),obj.S2(1:obj.inter_end));
            Sn((obj.inter_end+1):size(Sn,1)) = obj.S;%((obj.inter_end+1):size(Sn,1));
            H_mat = sparse(obj.I,obj.J,Sn,obj.mat_size,obj.mat_size);
   
        end
        
        
        
        
       
        %%
        function ind = center_index(obj,s_n,k)  % s_n is the sheet number
            % k is the basis index
            ind = 0;
            if s_n == 1
                ind = obj.index_inv(obj.sheet1.Origin_Index(k,1),1);
            elseif s_n == 2
                ind = obj.index_inv(obj.sheet2.Origin_Index(k,1),2);
            end
        end

        function v = center_vector(obj,s_n,k)
            v = zeros(obj.mat_size,1);
            v(obj.center_index(s_n,k))= 1;
        end
        
        
        
        function view(obj,color1,color2,c,col,shift,v_size)
            if ~exist('v_size','var')
                v_size = .02;
            end
            hold on
            %             I1 = obj.index(:,2) == 1;
            %             I2 = obj.index(:,2) == 2;
            %             I1 = obj.index(I1,1) ~= 0;
            %             I2 = obj.index(I2,2) ~= 0;
            
            
            s1_size = obj.sheet1.Basis_Size*obj.sheet1.Sheet_Size(1,1)*obj.sheet1.Sheet_Size(1,2);
            s2_size = obj.sheet2.Basis_Size*obj.sheet2.Sheet_Size(1,1)*obj.sheet2.Sheet_Size(1,2);
            temp_pos11 = obj.P(1:s1_size,1)' .*(obj.sheet1.Atom_Positions(1,:)-obj.sheet1.Origin(1,1))+shift(1,1);
            temp_pos12 = obj.P(1:s1_size,1)' .*(obj.sheet1.Atom_Positions(2,:)-obj.sheet1.Origin(2,1))+shift(2,1);
            temp_pos21 = obj.P((s1_size+1):(s1_size+s2_size),1)' .*(obj.sheet2.Atom_Positions(1,:)-obj.sheet2.Origin(1,1));
            temp_pos22 = obj.P((s1_size+1):(s1_size+s2_size),1)' .*(obj.sheet2.Atom_Positions(2,:)-obj.sheet2.Origin(2,1));
            scatter(temp_pos11,temp_pos12,v_size,color1);
            scatter(temp_pos21,temp_pos22,v_size,color2,'*');
            
            
            %scatter(obj.sheet1.Atom_Positions(1,obj.index(I1,1))-obj.sheet1.Origin(1,1),obj.sheet1.Atom_Positions(2,obj.index(I1,1))-obj.sheet1.Origin(2,1),color1);
            %scatter(obj.sheet2.Atom_Positions(1,obj.index(I2,1))-obj.sheet2.Origin(1,1),obj.sheet2.Atom_Positions(2,obj.index(I2,1))-obj.sheet2.Origin(2,1),color2);
            
            
            if c == 1
                for k = 1:obj.sheet1.Basis_Size
                    scatter(shift(1,1)+obj.sheet1.Atom_Positions(1,obj.sheet1.Origin_Index(k,1))-obj.sheet1.Origin(1,1),shift(2,1)+obj.sheet1.Atom_Positions(2,obj.sheet1.Origin_Index(k,1))-obj.sheet1.Origin(2,1),v_size,col);
                end
            elseif c == 2
                for k = 1:obj.sheet2.Basis_Size
                    scatter(obj.sheet2.Atom_Positions(1,obj.sheet2.Origin_Index(k,1))-obj.sheet2.Origin(1,1),obj.sheet2.Atom_Positions(2,obj.sheet2.Origin_Index(k,1))-obj.sheet2.Origin(2,1),v_size,col,'*');
                end
            end
        end

        
        function t = inter_func(obj,x,y,orbit1_i,orbit2_i,s1,s2) % o and m are orbitals
            
            orbit1 = (s1==1).*(orbit1_i)+(s1==2).*(orbit2_i);
            orbit2 = (s2==1).*(orbit2_i)+(s2==2).*(orbit1_i);
            orbit1 = (orbit1);
            orbit2 = (orbit2);
            
            x = (s1==1).*x -(s1==2).*x;
            y = (s1==1).*y -(s1==2).*y;
            
            r = sqrt(x.^2+y.^2+.00001);
            r_cut = 6;
            r_cut2 = 5;
            
            ac = acos(x./r);
            
            ac = (y<0) .* (2*pi-ac) + (y>=0).*ac;
            
            theta21 = (orbit1==1).*(ac + pi/6) + (orbit1==2).*(ac-pi/6);
            
            pi2o3 = 2*pi/3;
            theta21 = mod(theta21,pi2o3);
            
            ac = ac - obj.theta;
            
            theta12 = (orbit2==1).*(ac + pi/6) +(orbit2==2).*(ac-pi/6) ;
            
            
            theta12 = mod(theta12,pi2o3);
            
            % a = 2.46, lambda0 = .3155, lambda3 = -.0688, lambda6 = -.0083;
            % ci0 = 1.7543, ci3 = 3.4692, ci6 = 2.8764
            % x3 = .5212, x6 = 1.5206
            % k0 = 2.001, k6 = 1.5731
            
            V0 = .3155*exp(-1.7543*(r/2.46).^2).*cos(2.001*r./2.46);
            V3 = -.0688*(r/2.46).^2.*exp(-3.4692*(r/2.46 - .5212).^2);
            V6 = -.0083*exp(-2.8764*(r/2.46-1.5206).^2).*sin(1.5731*r/2.46);
            
            t = V0+V3.*(cos(3*theta12)+cos(3*theta21)) + V6.*(cos(6*theta12)+cos(6*theta21));
            
            t = (r > r_cut2).*t.*exp(1/(r_cut2-r_cut)^2-1./(r-r_cut).^2) + (r <= r_cut2).*t;
            % if isnan(t)
            %     'found a nan'
            %     [x,y,ot,orbit2_i,s1,s2]
            %     t=0;
            % end
            
        end
    end
end