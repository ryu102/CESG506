% Part 1
A = [1,1,2;1,4,1;1,4,4];
a = A\[3,4;5,1;6,3];

% Part 3
W = [3;-1];
F = [7/9, 3/9; -7/9, 6/9];
C = F'*F; 
[vectors,values] = eig(C);
N_eigenvalues = [values(1,1);values(2,2)];
N1 = vectors(:,1);
N2 = vectors(:,2);
U = 0;
for i = 1:2
    U = U + sqrt(N_eigenvalues(i))*(vectors(:,i)*vectors(:,i)');
end

R = F*inv(U);
V = F*inv(R);
[n,lambdan] = eig(V);
n_eigenvalues = [lambdan(1,1);lambdan(2,2)];
n1 = n(:,1);
n2 = n(:,2);
Check = R'*R;
w_right = R*U*W;
w_left = V*R*W;

%Part 4
I = eye(2,2);
E = .5*(F'*F-I);
e = .5*(I-(inv(F)'*inv(F)));
[c,g] = eig(E);
[f,h] = eig(e);
eig_E = [g(1,1);g(2,2)];
eig_e = [h(1,1);h(2,2)];


%Part 5
Test1 = (E - eig_E(1)*I)*N1;
Test2 = (E - eig_E(2)*I)*N2;
Test3 = (e - eig_e(1)*I)*n1;
Test4 = (e - eig_e(2)*I)*n2;

%Part 6
check_n1 = R*N1;
check_n2 = R*N2;

%Part 7
e_push = inv(F)'*E*inv(F);

%Part 8
E_pull = F'*e*F;