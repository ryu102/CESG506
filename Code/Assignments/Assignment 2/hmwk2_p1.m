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

%initial values
Pcr = [0;-.999];
ek = [0,1];
R = [0;0];
delta_uo = -inv(KTtot)*R;
u1 = -inv(KTtot)*Pcr;
u = [0;0];
g = ek*(u)+.01;
delta_lambda = -(g+dot(ek,delta_uo))/dot(ek,u1);
lambda = delta_lambda;

% storage
tol = .00001;
u_tot = zeros();
v_tot = zeros();
gamma = zeros();
Fx = zeros();
Fy = zeros();


for i = 1:120
    u(2) = -.01*i;
    v_tot(i) = u(2);
    
    %tried while loop but for loop was easier to implement
    for j = 1:100
        l1_vec = L1_vec + u';
        l2_vec = L2_vec + u';
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

        R = Ftot' + Pcr*lambda;

        Rtot = [R;g];
        
        %end loop if requirements are met
        if norm(Rtot) < tol
            u_tot(i) = u(1);
            gamma(i) = -lambda;
            Fx(i) = Ftot(1);
            Fy(i) = Ftot(2);
            break
        else
            
        I = eye(2);
        KT11 = EA/l1 * (n1'*n1);
        KT12 = f1/l1 * (I-n1'*n1);
        KT21 = EA/l2 *(n2'*n2);
        KT22 = f2/l2 * (I - n2'*n2);
        KTtot = KT11 + KT12 + KT21 + KT22;

        delta_uo = -inv(KTtot)*R;
        u1 = -inv(KTtot)*Pcr;
        g = ek*u + .01*i;
        delta_lambda = -(g+ek*delta_uo)/(ek*u1);
        lambda = lambda + delta_lambda;
        u = u + delta_uo + (u1*delta_lambda);

        k = norm(Rtot);
        end
    end

end

Fx_contour = [Fx',Fx'];
figure(1)
plot(v_tot,gamma); hold on
plot(u_tot,gamma)
title('Load Factor vs Displacement')
xlabel('Displacement (m)')
ylabel('Load Factor')
legend({'Vertical Displacement','Horizontal Displacement'})

figure(2)
plot(u_tot,v_tot)
title('Tracing loaded node')
xlabel('Horizontal Displacement(m)')
ylabel('Vertical Displacement(m)')







