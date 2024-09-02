function t = inter_graphene(x,y,orbit1,orbit2,theta)
            
            % orbit1 = (s1==1).*orbit1_i+(s1==2).*orbit2_i;
            % orbit2 = (s2==1).*orbit2_i+(s2==2).*orbit1_i;
            % 
            % x = (s1==1).*x -(s1==2).*x;
            % y = (s1==1).*y -(s1==2).*y;

            r = sqrt(x.^2+y.^2+.00001);
            r_cut = 6;
            r_cut2 = 5;
            
            ac = acos(x./r);
            
            ac = (y < 0) .* (2*pi-ac) + ( y >= 0).*ac;
            
            theta21 = (orbit1==1).*(ac + pi/6) + (orbit1==2).*(ac-pi/6);
            
            pi2o3 = 2*pi/3;
            theta21 = mod(theta21,pi2o3);
            
            ac = ac - theta;
            
            theta12 = (orbit2==1).*(ac + pi/2) +(orbit2==2).*(ac-pi/6) ;
            
            
            %theta12 = mod(theta12,pi2o3);
            
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