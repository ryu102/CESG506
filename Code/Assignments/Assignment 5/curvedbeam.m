function [Fe, Ke] = curvedbeam(Xi,Xj, ui, uj, EA, EI)
    
    L = Xj(1) - Xi(1);
    hi = Xi(2);
    hj = Xj(2);
    dh = (hj - hi)/L;

    qu = [ui(1);uj(1)];
    qv = [ui(2); ui(3); uj(2); uj(3)];

    Fe = zeros(6,1);
    Ke = zeros(6,6);

    % gauss 
    gauss_x = [0.211, 0.788];
    weight = [.5, .5];

    for i = 1:2

        x =gauss_x(i);
        W = weight(i)*L;

        dNu = [-1,1]/L;
        dNv = [-6*x + 6*x^2, L*(1-4*x+3*x^2),...
              (6*x-6*x^2), L*(-2*x+3*x^2)]/L;
        ddNv = [-6 + 12*x, L*(-4+6*x),...
               (6-12*x), L*(-2+6*x)]/(L*L);

        du0 = dNu * qu;
        dv0 = dNv * qv;
        phi = ddNv * qv;
        eps0 = du0 + dh*dv0 + 0.5*dv0^2;

        f = EA*eps0;
        M = EI*phi;

        B_phi = [0, ddNv(1), ddNv(2), 0, ddNv(3), ddNv(4)];
        B_eps = [dNu(1), (dh+dv0)*dNv(1), (dh+dv0)*dNv(2),...
                dNu(2), (dh+dv0)*dNv(3), (dh+dv0)*dNv(4)];

        dNv_mat = [0,dNv(1), dNv(2),0, dNv(3), dNv(4)];

        Fe = Fe + (f*B_eps' + M*B_phi')*W;
        Ke = Ke + (EA*B_eps'*B_eps + f*dNv_mat'*dNv_mat + EI*B_phi'*B_phi)*W;
    end 
end
    
    