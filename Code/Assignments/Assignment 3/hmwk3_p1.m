EA = 2100;
L_vec = [5.5, 0.5];
L = sqrt(dot(L_vec,L_vec));
n = L_vec/L;
K = EA/L * L_vec(2)*L_vec(2);
Pcr = -.3;

R = 0;
vo = K\R;
v1 = K\Pcr;
vn = 0;
gamman = 0;
g = 0;
alpha = 0.1;
s=0;
iterations = 1;

v_tot = zeros();
gamma_tot = zeros();
F = zeros();
arclength = zeros();
normR = zeros();
Residuals = zeros();
normRspec = zeros();

gamma_tot(1) = gamman;
v_tot(1) = vn;


for i = 2:35
    if i == 2
        gamma = gamman + .25;
        v = vn + v1*.25;
        ds2 = (v-vn)^2 + alpha*(gamma-gamman)^2;
        ds = sqrt(ds2);
    else
        v = 2*v_tot(end) - v_tot(end-1);
        gamma = 2*gamma_tot(end) - gamma_tot(end-1);
    end
    
    for j = 1:100
        
        l_vec = L_vec + [0,v];
        l = sqrt(dot(l_vec,l_vec));
        n = l_vec/l;
        lambda = l/L;
        strain = log(lambda);
        f = EA*strain;
        Ftot = f*n(2);

        R = -Ftot + Pcr*gamma;
        
        change_v = v-vn;
        change_gamma = gamma - gamman;
        g = -ds2 + change_v^2 + alpha*change_gamma^2; 
    
        Rtot = [R;g];
        check = norm(Rtot);
        Residuals(j,i) = check;
        
        iterations = iterations + 1;
        normRspec(iterations) = check;
        

        if norm(Rtot) < .0000001
           break
        else
        end
        k1 = EA/l * (n(2)*n(2));
        k2 = f/l * (1-n(2)*n(2));
        ktot = k1+ k2;

        vo = ktot\R;
        v1 = ktot\Pcr;
        dgamma = -(g+2*change_v*vo)/(2*change_v*v1 + 2*alpha*change_gamma);
        dv = vo + v1*dgamma;
        gamma = gamma + dgamma;
        v = v + dv;
        
       
    end
    v_tot(i) = v;
    gamma_tot(i) = gamma;
    F(i) = Ftot;
    arclength(i) = s+ds;
    normR(i) = check;
    s = s + ds;
    gamman = gamma;
    vn = v;
    K = ktot;
end

% done in excel
uhw1 = linspace(0,-1.2,13);
fhw1 = [0,-.225,-.30,-.263,-.15,0,.150,.263,.301,.225,0,-.409,-1.03];
pltfhw1 =-fhw1/.3;

iter = linspace(0,iterations,iterations);

figure(1)
plot(v_tot,gamma_tot); hold on
scatter(uhw1,pltfhw1)
title('Load Factor vs Displacement')
xlabel('Displacement (m)')
ylabel('Load Factor')
legend({'Hmwk 3','Hmwk 1'})


figure(2)
semilogy(iter,normRspec)
title('|R| vs iterations')
xlabel('iterations')
ylabel('|R|')

figure(3)
semilogy(arclength,normR)
title('|R| vs arclength')
xlabel('arclength')
ylabel('|R|')




        