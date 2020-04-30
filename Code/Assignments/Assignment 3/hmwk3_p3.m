h = 5/11; % story height
p = 6.2/2; % calculated load
nodes = [1  0   0;
         2  .25 0;
         3  0   h;
         4  .25 h;
         5  0   2*h;
         6  .25 2*h;
         7  0   3*h;
         8  .25 3*h;
         9  0   4*h;
         10 .25 4*h;
         11 0   5*h;
         12 .25 5*h;
         13 0   6*h;
         14 .25 6*h;
         15 0   7*h;
         16 .25 7*h;
         17 0   8*h;
         18 .25 8*h;
         19 0   9*h;
         20 .25 9*h;
         21 0   10*h;
         22 .25 10*h;
         23 0   11*h;
         24 .25 11*h;];
     
mesh = [1   1   3   2000;
        2   1   4   5000;
        3   2   4   2000;
        4   3   4   5000;
        5   3   5   2000;
        6   3   6   5000;
        7   4   6   2000;
        8   5   6   5000;
        9   5   7   2000;
        10  5   8   5000;
        11  6   8   2000;
        12  7   8   5000;
        13  7   9   2000;
        14  7   10  5000;
        15  8   10  2000;
        16  9   10  5000;
        17  9   11  2000;
        18  9   12  5000;
        19  10  12  2000;
        20  11  12  5000;
        21  11  13  2000;
        22  11  14  5000;
        23  12  14  2000;
        24  13  14  5000;
        25  13  15  2000;
        26  13  16  5000;
        27  14  16  2000;
        28  15  16  5000;
        29  15  17  2000;
        30  15  18  5000;
        31  16  18  2000;
        32  17  18  5000;
        33  17  19  2000;
        34  17  20  5000;
        35  18  20  2000;
        36  19  20  5000;
        37  19  21  2000;
        38  19  22  5000;
        39  20  22  2000;
        40  21  22  5000;
        41  21  23  2000;
        42  21  24  5000;
        43  22  24  2000;
        44  23  24  5000;];
    
support = [1    1   1;
           2    1   1;];

P = [0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   -p;
     0   -p;]';
 
U = [0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;]';
 
 idxFree = linspace(3,48,46)';
 idxFixed= [1;2;3;4];
 
[Pcr, Ft, Kt] = Assemble(nodes, mesh, support, P, U,idxFree,idxFixed);
R = zeros(44,1);
vo = Kt\R;
v1 = Kt\Pcr;
vn = zeros(44,1);
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
   
for i = 2:116
    if i == 2
        gamma = gamman + .5;
        v = vn + v1*.5;
        ds2 = dot(v-vn,v-vn) + alpha*(gamma-gamman)^2;
        ds = sqrt(ds2);
    else
        v = 2*v_tot{end} - v_tot{end-1};
        gamma = 2*gamma_tot(end) - gamma_tot(end-1);
    end
    
    for j = 1:10

        Ux = [0;0;v(1:2:end)];
        Uy = [0;0;v(2:2:end)];
        U = [Ux,Uy]';

        [Pe, Ftot, Ktot] = Assemble(nodes, mesh, support, P, U,idxFree,idxFixed);

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
plot(v_values(42,:),gamma_tot); hold on
plot(v_values(41,:),gamma_tot)
title('Top Left Node Displacement vs Load Factor')
xlabel('Displacement (m)')
ylabel('Load Factor')
legend({'Vertical Displacement','Horizontal Displacement'})

figure(2)
plot(v_values(44,:),gamma_tot); hold on
plot(v_values(43,:),gamma_tot)
title('Top Right Node Displacement vs Load Factor')
xlabel('Displacement (m)')
ylabel('Load Factor')
legend({'Vertical Displacement','Horizontal Displacement'})

% to make plot start from 5m at top
top_dis = ones(1,i)*5;
disp_left = top_dis + v_values(44,:);
disp_right = top_dis + v_values(42,:);

figure(3)
plot(v_values(43,:),disp_left); hold on
plot(v_values(41,:),disp_right)
title('Top Node Paths')
xlabel('Horizontal Displacement (m)')
ylabel('Vertical Displacement (m)')
legend({'Top Left Node','Top Right Node'})
grid on

