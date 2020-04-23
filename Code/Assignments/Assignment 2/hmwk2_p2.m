%inputs
EA =2100;

%establish nodes
node1 = [0,0,-5];
node2 = [9.5,0,-5];
node3 = [0,0,0];
node4 = [9.5,0,0];
node5 = [5.5,0.5,-3.75];
node6 = [5.5,0.5,-1.25];

% establish member vector
mem16_vec = node6-node1;
mem15_vec = node5-node1;
mem25_vec = node5-node2;
mem56_vec = node5-node6;
mem45_vec = node5-node4;
mem46_vec = node6-node4;
mem36_vec = node6-node3;

% establish initial length
L_mem16 = sqrt(dot(mem16_vec,mem16_vec));
L_mem15 = sqrt(dot(mem15_vec,mem15_vec));
L_mem25 = sqrt(dot(mem25_vec,mem25_vec));
L_mem56 = sqrt(dot(mem56_vec,mem56_vec));
L_mem45 = sqrt(dot(mem45_vec,mem45_vec));
L_mem46 = sqrt(dot(mem46_vec,mem46_vec));
L_mem36 = sqrt(dot(mem36_vec,mem36_vec));

%establish initial normal vectors
N_mem16 = mem16_vec/L_mem16;
N_mem15 = mem15_vec/L_mem15;
N_mem25 = mem25_vec/L_mem25;
N_mem56 = mem56_vec/L_mem56;
N_mem45 = mem45_vec/L_mem45;
N_mem46 = mem46_vec/L_mem46;
N_mem36 = mem36_vec/L_mem36;

%establish inital K matrix
K_mem16 = EA/L_mem16 * (N_mem16'*N_mem16);
K_mem15 = EA/L_mem15 * (N_mem15'*N_mem15);
K_mem25 = EA/L_mem25 * (N_mem25'*N_mem25);
K_mem56 = EA/L_mem56 * (N_mem56'*N_mem56);
K_mem45 = EA/L_mem45 * (N_mem45'*N_mem45);
K_mem46= EA/L_mem46 * (N_mem46'*N_mem46);
K_mem36 = EA/L_mem36 * (N_mem36'*N_mem36);
K11 = K_mem15+K_mem25+K_mem45+K_mem56;
K12 = -K_mem56;
K21 = -K_mem56;
K22 = K_mem16+K_mem36+K_mem46+K_mem56;
Ktot = [K11, K12; K12, K22];

%initial values
Pcr = [0;-.999;0;0;0;0];
ek = [0,1,0,0,0,0];
R = [0;0;0;0;0;0];
delta_uo = -inv(Ktot)*R;
u1 = -inv(Ktot)*Pcr;
u = [0;0;0;0;0;0];
g = ek*u+.01;
delta_lambda = -(g*dot(ek,delta_uo))/dot(ek,u1);
lambda = delta_lambda;

%storage
tol = .00001;
u5 = zeros();
v5 = zeros();
z5 = zeros();
u6 = zeros();
v6 = zeros();
z6 = zeros();
gamma = zeros();

planeu5 = zeros();
planev5 = zeros();
planez5 = zeros();
planeu6 = zeros();
planev6 = zeros();
planez6 = zeros();

Fx = zeros();
Fy = zeros();
Fz = zeros();

for i = 1:120
    u(2) = -.01*i;
    v5(i) = u(2);
    planev5(i) = node5(2) - u(2);
    
    for j = 1:10000
        lmem15_vec = mem15_vec + u(1:3)';
        lmem16_vec = mem16_vec + u(4:6)';
        lmem25_vec = mem25_vec + u(1:3)';
        lmem36_vec = mem36_vec + u(4:6)';
        lmem45_vec = mem45_vec + u(1:3)';
        lmem46_vec = mem46_vec + u(4:6)';
        lmem56_vec = mem56_vec + u(1:3)'-u(4:6)';
        
        l_mem15 = sqrt(dot(lmem15_vec,lmem15_vec));
        l_mem16 = sqrt(dot(lmem16_vec,lmem16_vec));
        l_mem25 = sqrt(dot(lmem25_vec,lmem25_vec));
        l_mem36 = sqrt(dot(lmem36_vec,lmem36_vec));
        l_mem45 = sqrt(dot(lmem45_vec,lmem45_vec));
        l_mem46 = sqrt(dot(lmem46_vec,lmem46_vec));
        l_mem56 = sqrt(dot(lmem56_vec,lmem56_vec));
        
        n_mem15 = lmem15_vec/l_mem15;
        n_mem16 = lmem16_vec/l_mem16;
        n_mem25 = lmem25_vec/l_mem25;
        n_mem36 = lmem36_vec/l_mem36;
        n_mem45 = lmem45_vec/l_mem45;
        n_mem46 = lmem46_vec/l_mem46;
        n_mem56 = lmem56_vec/l_mem56;
        
        lambda15 = l_mem15/L_mem15;
        lambda16 = l_mem16/L_mem16;
        lambda25 = l_mem25/L_mem25;
        lambda36 = l_mem36/L_mem36;
        lambda45 = l_mem45/L_mem45;
        lambda46 = l_mem46/L_mem46;
        lambda56 = l_mem56/L_mem56;
        
        strain15 = log(lambda15);
        strain16 = log(lambda16);
        strain25 = log(lambda25);
        strain36 = log(lambda36);
        strain45 = log(lambda45);
        strain46 = log(lambda46);
        strain56 = log(lambda56);
        
        f15 = EA*strain15;
        f16 = EA*strain16;
        f25 = EA*strain25;
        f36 = EA*strain36;
        f45 = EA*strain45;
        f46 = EA*strain46;
        f56 = EA*strain56;
        
        F5_tot = f15*n_mem15 + f25*n_mem25 + f45*n_mem45 + f56*n_mem56;
        F6_tot = f16*n_mem16 + f36*n_mem36 + f46*n_mem46 - f56*n_mem56;
        
        R5 = lambda*Pcr(1:3) + F5_tot';
        R6 = lambda*Pcr(4:6) + F6_tot'; 
        R = [R5;R6];
        Rtot = [R;g];
        check = norm(Rtot);
        
        if norm(Rtot) < tol
            u5(i) = u(1);
            planeu5(i) = node5(1) - u(1);
            z5(i) = u(3);
            planez5(i) = node5(3) - u(3);
            u6(i) = u(4);
            planeu6(i) = node6(1) - u(4);
            v6(i) = u(5);
            planev6(i) = node6(2) - u(5);
            z6(i) = u(6);
            planez6(i) = node6(3) - u(6);
            gamma(i) = -lambda;
            break
        else
        
        % create delta values to add to existing disp and lambda
        I = eye(3);
        k_mem16 = EA/l_mem16 * (n_mem16'*n_mem16)+ f16/l_mem16*(I-n_mem16'*n_mem16);
        k_mem15 = EA/l_mem15 * (n_mem15'*n_mem15)+ f15/l_mem15*(I-n_mem15'*n_mem15);
        k_mem25 = EA/l_mem25 * (n_mem25'*n_mem25)+ f25/l_mem25*(I-n_mem25'*n_mem25);
        k_mem56 = EA/l_mem56 * (n_mem56'*n_mem56)+ f56/l_mem56*(I-n_mem56'*n_mem56);
        k_mem45 = EA/l_mem45 * (n_mem45'*n_mem45)+ f45/l_mem45*(I-n_mem45'*n_mem45);
        k_mem46 = EA/l_mem46 * (n_mem46'*n_mem46)+ f46/l_mem46*(I-n_mem46'*n_mem46);
        k_mem36 = EA/l_mem36 * (n_mem36'*n_mem36)+ f36/l_mem36*(I-n_mem36'*n_mem36);
        
        k11 = k_mem15+k_mem25+k_mem45+k_mem56;
        k12 = -k_mem56;
        k22 = k_mem16+k_mem36+k_mem46+k_mem56;
        ktot = [k11, k12; k12, k22];
        
        delta_uo = -inv(ktot)*R;
        u1 = -inv(ktot)*Pcr;
        g = ek*u + .01*i;
        delta_lambda = -(g+ek*delta_uo)/(ek*u1);
        lambda = delta_lambda + lambda;
        u = u + delta_uo + (u1*delta_lambda);
        end
    end
end

%node 6
figure(1)
plot(v5,gamma); hold on
plot(u5,gamma); hold on
plot(z5,gamma); hold on
plot(v6,gamma); hold on
plot(u6,gamma); hold on
plot(z6,gamma)
title('Load Factor vs Displacement')
xlabel('Displacement (m)')
ylabel('Load Factor')
legend({'Node 5: Vertical Displacement','Node 5: Horizontal Displacement','Node 5: Out of Plane Displacement','Node 6: Vertical Displacement','Node 6: Horizontal Displacement','Node 6: Out of Plane Displacement'})

figure(2)
plot(u5,v5)
title('Tracing node 5 X&Y')
xlabel('Horizontal Displacement (m)')
ylabel('Vertical Displacement (m)')

figure(3)
plot(z5,v5)
title('Tracing node 5 X&Z')
xlabel('Out of Plane Disaplacement (m)')
ylabel('Vertical Displacement (m)')

figure(5)
plot(u6,v6)
title('Tracing node 6 X&Y')
xlabel('Horizontal Displacement (m)')
ylabel('Vertical Displacement (m)')

figure(6)
plot(z6,v6)
title('Tracing node 6 X&Z')
xlabel('Out of Plane Disaplacement(m)')
ylabel('Vertical Displacement (m)')

figure(7)
plot(planez5,planev5); hold on
plot(planez6,planev6)
title('Planar View Vertical')
xlabel('Out of Plane Location(m)')
ylabel('Vertical Location (m)')
legend({'Node 5', 'Node 6'})

figure(8)
plot(planez5,planeu5); hold on
plot(planez6,planeu6)
title('Planar View Horizontal')
xlabel('Out of Plane Location(m)')
ylabel('Horizontal Location (m)')
legend({'Node 5', 'Node 6'})
        
        
    


