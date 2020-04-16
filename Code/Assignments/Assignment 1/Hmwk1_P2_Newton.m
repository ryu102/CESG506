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
KT_L1 = EA/L1 *(N1'*N1);
KT_L2 = EA/L2 *(N2'*N2);
KTtot = KT_L1 + KT_L2;

Pcr = [0;.9786];
lambda = [0, 0.25, 0.5, 0.75, 0.99, 0.999];

%Initial Values and storage
Residual = zeros(10,5);
u_change = [0,0];
Ftot = [0,0];
tol = .000000001;

for i = 1:5  
    timestep = 0;
    R = lambda(i+1)*Pcr + Ftot';
    h = -inv(KTtot)*R;
    Residual(timestep+1,i) = norm(R);
    
    while norm(R) > tol
        timestep = timestep +1;
        u_change = u_change + h';
        
        l1_vec = L1_vec + u_change;
        l2_vec = L2_vec + u_change;
        l1 = sqrt(dot(l1_vec,l1_vec));
        l2 = sqrt(dot(l2_vec,l2_vec));
        n1 = l1_vec/l1;
        n2 = l2_vec/l2;
        lambda1 = l1/L1;
        lambda2 = l2/L2;
        strain1 = log(lambda1);
        strain2 = log(lambda2);
        f1 = EA*strain1;
        f2 = EA*strain2;
        Ftot = f1*n1 + f2*n2;
        
        R = lambda(i+1)*Pcr + Ftot';
        Residual(timestep+1,i) = norm(R);
        
        I = eye(2);
        KT11 = EA/l1 * (n1'*n1);
        KT12 = f1/l1 * (I-n1'*n1);
        KT21 = EA/l2 *(n2'*n2);
        KT22 = f2/l2 * (I - n2'*n2);
        KTtot = KT11 + KT12 + KT21 + KT22;
        
        h = -inv(KTtot)*R;
    end
end





