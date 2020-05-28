% set initial conditions
EA = 10000;
EI = 10;
num_nodes = 5;
U = zeros(3,num_nodes);
L = 5;
idxFixed = [1,2,3*num_nodes-2,3*num_nodes-1];
idxFree = [3:1:3*(num_nodes-1),3*num_nodes]; 

nodes= zeros(num_nodes,3);
mesh = zeros(num_nodes-1,5);
P = zeros(3,num_nodes);
nodes(1,1) =1;

% create nodes matrix
for i = 2:num_nodes
    nodes(i,1) = i;
    X = L/(num_nodes-1)*(i-1);
    nodes(i,2) = X;
    nodes(i,3) = -.08*X^2 + .4*X;   
end

% create mesh matrix
for i = 1:(num_nodes-1)
    mesh(i,1) = i;
    mesh(i,2) = i;
    mesh(i,3) = i+1;
    mesh(i,4) = EA;
    mesh(i,5) = EI;
end

%force vector
for i = 1:num_nodes-1
    x1 = nodes(i,2);
    x2 = nodes(i+1,2);
    P_node1 = .004*x1+.99;
    P_node2 = .004*x2+.99;
    P_tot = P_node1*(x2-x1) + .5*(x2-x1)*(P_node2-P_node1);
    P(2,i) = P(2,i)-.5*P_tot;
    P(2,i+1) = P(2,i+1)-.5*P_tot;
end

[Pcr, Ft, Kt] = Assemble(nodes, mesh, P, U,idxFree,idxFixed);

R = zeros(num_nodes*3-4,1);
vo = Kt\R;
v1 = Kt\Pcr;
vn = zeros(num_nodes*3-4,1);
gamman = 0;
g = 0;
alpha = 0.1;
s=0;

v_tot = {};
gamma_tot = zeros();
Fy = zeros();
arclength = zeros();
normR = zeros();
Residuals = zeros();

gamma_tot(1) = gamman;
v_tot{1} = vn;
   
%arc length
for i = 2:39
    if i == 2
        gamma = gamman + .2;
        v = vn + v1*.2;
        ds2 = dot(v-vn,v-vn) + alpha*(gamma-gamman)^2;
        ds = sqrt(ds2);
    else
        v = 2*v_tot{end} - v_tot{end-1};
        gamma = 2*gamma_tot(end) - gamma_tot(end-1);
    end
    
    for j = 1:10

        Ux = [0;v(2:3:end-1);0];
        Uy = [0;v(3:3:end-1);0];
        theta = [v(1:3:end);v(end)];
        U = [Ux,Uy,theta]';

        [Pe, Ftot, Ktot] = Assemble(nodes, mesh, P, U,idxFree,idxFixed);

        R = -Ftot + Pcr*gamma;

        change_v = v-vn;
        change_gamma = gamma - gamman;
        g = -ds2 + dot(change_v,change_v) + alpha*change_gamma^2;

        Rtot = [R;g];
        check = norm(Rtot);
        Residuals(j,i) = check;

        if norm(Rtot) < 1e-10
           break
        end

        vo = Ktot\R;
        v1 = Ktot\Pcr;
        dgamma = -(g+2*dot(change_v,vo))/(2*dot(change_v,v1) + 2*alpha*change_gamma);
        dv = vo + v1*dgamma;
        gamma = gamma + dgamma;
        v = v + dv;

    end
    v_tot{i} = v;
    gamma_tot(i) = gamma;
    arclength(i) = s+ds;
    normR(i) = check;
    s = s + ds;
    gamman = gamma;
    vn = v;
end

v_values = cell2mat(v_tot);

figure(1)
plot(v_values(2,:),gamma_tot); hold on
plot(v_values(3,:),gamma_tot); hold on
plot(v_values(5,:),gamma_tot); hold on
plot(v_values(6,:),gamma_tot); hold on
plot(v_values(8,:),gamma_tot); hold on
plot(v_values(9,:),gamma_tot); hold on
title('Load Factor vs Displacement')
xlabel('Displacement (m)')
ylabel('Load Factor')
legend({'Horizontal Displacement at 1/4 PT','Vertical Displacement at 1/4 PT',...
        'Horizontal Displacement at 1/2 PT','Vertical Displacement at 1/2 PT',...
        'Horizontal Displacement at 3/4 PT','Vertical Displacement at 3/4 PT'},'location','northwest')

Ux_final = nodes(:,2)+U(1,:)';
Uy_final = nodes(:,3)+U(2,:)';

figure(2)
plot(nodes(:,2),nodes(:,3)); hold on
plot(Ux_final,Uy_final)
title('Undeformed vs Deformed Shape')
xlabel('Horizontal Distance (m)')
ylabel('Vertical Distance (m)')
legend({'Undeformed Shape','Deformed Shape'})
grid on


