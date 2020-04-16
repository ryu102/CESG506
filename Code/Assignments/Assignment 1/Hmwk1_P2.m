% Inputs 
EA = 2100;
W1 = 5.5;
W2 = -4;
H = 0.5;
L1_vec = [W1,H];
L2_vec = [W2,H];
L1 = sqrt(dot(L1_vec,L1_vec));
L2 = sqrt(dot(L2_vec,L2_vec));
N1 = L1_vec/L1;
N2 = L2_vec/L2;

%Solve in x direction for position vectors
V = linspace(0,2*H,100);
u_actual = zeros(1,100);
Fx = zeros(100,250);
Fx(:,1) = 1;

for i = 1:100
    for j = 2:250
        ux = j/5000;
        ux_change = [-ux,-V(i)];
        l1_vecx = L1_vec + ux_change;
        l2_vecx = L2_vec + ux_change;
        l1x = sqrt(dot(l1_vecx,l1_vecx));
        l2x = sqrt(dot(l2_vecx,l2_vecx));
        n1x = l1_vecx/l1x;
        n2x = l2_vecx/l2x;
        lambda1x = l1x/L1;
        lambda2x = l2x/L2;
        strain1x = log(lambda1x);
        strain2x = log(lambda2x);
        f1x = EA*strain1x;
        f2x = EA*strain2x;
        F = f1x*n1x + f2x*n2x;
        Fx(i,j) = F(:,1);
        
        if Fx(i,j)*Fx(i,j-1) < 0 
            u_actual(i) = -ux;
        else
        end
    end
end

% Force in Y to get Forces
Fy = zeros(1,100);
for k = 1:100
    uy_change = [u_actual(k),-V(k)];
    l1_vecy = L1_vec + uy_change;
    l2_vecy = L2_vec + uy_change;
    l1y = sqrt(dot(l1_vecy,l1_vecy));
    l2y = sqrt(dot(l2_vecy,l2_vecy));
    n1y = l1_vecy/l1y;
    n2y = l2_vecy/l2y;
    lambda1y = l1y/L1;
    lambda2y = l2y/L2;
    strain1y = log(lambda1y);
    strain2y = log(lambda2y);
    f1y = EA*strain1y;
    f2y = EA*strain2y;
    f_tot = f1y*n1y + f2y*n2y;
    Fy(k) = f_tot(:,2);
end

plot(V,Fy)
xlabel("Vertical Displacement (m)")
ylabel("Force (KN)")
Pcr = max(Fy);
    
    
    
    
        
        