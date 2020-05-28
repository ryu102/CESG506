% set initial conditions
EA = 10000;
EI = 10;
num_nodes = 17;
U = zeros(3,num_nodes);
L = 5;
idxFixed = [1,2,3*num_nodes-1];
idxFree = [3:1:3*(num_nodes)-2,3*num_nodes]; 

nodes= zeros(num_nodes,3);
mesh = zeros(num_nodes-1,5);
P = zeros(3,num_nodes);
P(1,num_nodes) = -3.95;
nodes(1,1) =1;

% create nodes matrix
for i = 2:num_nodes
    nodes(i,1) = i;
    X = L/(num_nodes-1)*(i-1);
    nodes(i,2) = X;
    nodes(i,3) = 0;   
end

% create mesh matrix
for i = 1:(num_nodes-1)
    mesh(i,1) = i;
    mesh(i,2) = i;
    mesh(i,3) = i+1;
    mesh(i,4) = EA;
    mesh(i,5) = EI;
end

[Pcr, Ft, Kt] = Assemble(nodes, mesh, P, U,idxFree,idxFixed);

R = zeros(num_nodes*3-3,1);
vo = Kt\R;
v1 = Kt\Pcr;
vn = zeros(num_nodes*3-3,1);
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
determinants = zeros();

gamma_tot(1) = gamman;
v_tot{1} = vn;
   
%arc length
for i = 2:65
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

        Ux = [0;v(2:3:end)];
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
        determinant = min(eig(Ktot));

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
    determinants(i) = determinant;
    s = s + ds;
    gamman = gamma;
    vn = v;
    
    if determinant < 0
        break
    end
end

v_values = cell2mat(v_tot);

figure(1)
plot(gamma_tot,determinants)
title('Determinant of K vs Load Factor')
xlabel('Load Factor')
ylabel('Determinant of K')



