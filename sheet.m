classdef sheet
    properties
        Lattice
        Sheet_Size           % max sheet sizes (1,2) vector
        Basis_Size
        Basis_Positions
        Node_Positions
        Atom_Positions
        Origin               % this is considered the origin when taking input.
        index                % (i,j,b) coordinates to natural numbers coordinates
        index_inv            % opposite of above
        Ratio                % used in neighbor calculations
        Origin_Index
        
    end
    methods
        function obj = sheet(L,Sheet_Size,basis_positions,shift) % constructor
            if ~exist('shift','var')
                shift = [0;0];
            end
            
            obj.Lattice = L;
            obj.Sheet_Size = Sheet_Size;
            obj.Basis_Positions = basis_positions;
            obj.Basis_Size = size(basis_positions,2);
            obj.Node_Positions = zeros(2,Sheet_Size(1,1)*Sheet_Size(1,2));
            obj.Atom_Positions = zeros(2,Sheet_Size(1,1)*Sheet_Size(1,2)*obj.Basis_Size);
            obj.Origin_Index = zeros(obj.Basis_Size,1);
            %s_size = Sheet_Size(1,1)*Sheet_Size(1,2)*obj.Basis_Size;
            
            obj.Origin = L*(Sheet_Size'-1)*.5 + shift;
            obj.index = zeros(Sheet_Size(1,1),Sheet_Size(1,2),2);
            obj.index_inv = zeros(Sheet_Size(1,1)*Sheet_Size(1,2)*2,3);
            obj.Ratio = abs(dot(obj.Lattice(:,1),obj.Lattice(:,2)))/(norm(obj.Lattice(:,1))*norm(obj.Lattice(:,2)));
            obj.Ratio = min(obj.Ratio,1-obj.Ratio)^(-1);
            for i = 1:Sheet_Size(1,1)           % build sheet
                for j = 1:Sheet_Size(1,2)
                    for k = 1:obj.Basis_Size
                        obj.index(i,j,k) = (i-1)*Sheet_Size(1,2)*obj.Basis_Size + (j-1)*obj.Basis_Size + k;
                        obj.index_inv(obj.index(i,j,k),1) = i;
                        obj.index_inv(obj.index(i,j,k),2) = j;
                        obj.index_inv(obj.index(i,j,k),3) = k;
                        
%                         if i == (Sheet_Size(1,1)-mod(Sheet_Size(1,1),2))/2 && j == (Sheet_Size(1,2)-mod(Sheet_Size(1,2),2))/2
%                             obj.Origin_Index(k,1) = obj.index(i,j,k);
%                         end
                        
                        obj.Node_Positions(:,(i-1)*Sheet_Size(1,2) + j) = L*[i-1;j-1];
                        obj.Atom_Positions(:,obj.index(i,j,k)) = L*[i-1;j-1]+basis_positions(:,k);
                        if norm(L*[i-1;j-1]-obj.Origin) < .5*max(norm(obj.Lattice(:,1)),norm(obj.Lattice(:,2)))
                            obj.Origin_Index(k,1) = obj.index(i,j,k);
                        end
                    end
                end
            end

        end

        function Neighbors = neighbors(obj,Position, Cutoff,P)
            noP = 0; % we assume there is a P
            if ~exist('P','var')
                noP = 1;
                P = ones(size(obj.Atom_Positions,2),1);
            end
            
            Cutoff_Adj = Cutoff*obj.Ratio;
            Cutoff_Squared = Cutoff^2;
            Position_Adj = Position + obj.Origin;
            coord = floor(obj.Lattice^(-1)*(Position_Adj));
            
            if coord(1,1) <= 0
                coord(1,1) = 1;
            end
            if coord(2,1) <= 0
                coord(2,1) = 1;
            end
            
            % not taking into account basis!
            
            min1 = floor(max(1,coord(1,1)-Cutoff_Adj*norm(obj.Lattice(:,1))^(-1)-1));
            max1 = ceil(min(obj.Sheet_Size(1,1),coord(1,1)+Cutoff_Adj*norm(obj.Lattice(:,1))^(-1)+1));
            min2 = floor(max(1,coord(2,1)-Cutoff_Adj*norm(obj.Lattice(:,2))^(-1)-1));
            max2 = ceil(min(obj.Sheet_Size(1,2),coord(2,1)+Cutoff_Adj*norm(obj.Lattice(:,2))^(-1)+1));
            % count = 0;
            %s_size = size(obj.Atom_Positions,2);
            %Position_Squared = zeros(s_size,1);
            


            % for j = 1:s_size
            %     Position_Squared(j,1) = (obj.Atom_Positions(1,j)-Position_Adj(1,1))^2 ...
            %         + (obj.Atom_Positions(2,j)-Position_Adj(2,1))^2;
            % end

            % [i,j,k] = meshgrid(min1:max1,min2:max2, 1:obj.Basis_Size);
            % i = permute(i,[2 1 3]);
            % j = permute(j,[2 1 3]);

            l = obj.index(min1:max1,min2:max2, 1:obj.Basis_Size);
            l = l(:);

            Position_Squared = (obj.Atom_Positions(1,l)-Position_Adj(1,1)).^2 ...
            + (obj.Atom_Positions(2,l)-Position_Adj(2,1)).^2;

            Neighbors = l(Position_Squared(:) < Cutoff_Squared & P(l)==1);

            % if Position_Squared(l,1)<Cutoff_Squared
            %     if noP == 0
            %         if P(l,1) == 0
            %             count = count - 1;
            %         end
            %     end
            %     count = count + 1;
            % end

            % Neighbors = zeros(count,1);
            % count_pos = 1;
            % for i = min1:max1
            %     for j = min2:max2
            %         for k = 1:obj.Basis_Size
            %             l = obj.index(i,j,k);
            % 
            %             if Position_Squared(l,1)<Cutoff_Squared
            %                 if noP == 1
            %                     Neighbors(count_pos,1) = l;
            %                     count_pos = count_pos+1;
            %                 end
            %                 if noP == 0
            %                 	if P(l,1) == 1
            %                         Neighbors(count_pos,1) = l;
            %                         count_pos = count_pos + 1;
            %                     end
            %                end
            %             end
            %         end
            %     end
            %end 
        end
        
        
        function view(obj,color,extra_command,color2)
            if ~exist('extra_command','var')
                extra_command = 'nothing';
            end

            scatter(obj.Atom_Positions(1,:)-obj.Origin(1,1),obj.Atom_Positions(2,:)-obj.Origin(2,1),color);
            if strcmp(extra_command,'highlight_center')
                hold on
                for k = 1:obj.Basis_Size
                    scatter(obj.Atom_Positions(1,obj.Origin_Index(k,1))-obj.Origin(1,1),obj.Atom_Positions(2,obj.Origin_Index(k,1))-obj.Origin(2,1),color2);
                end
            end
        end
        

    end
end
        