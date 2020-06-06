% Part 3
J = [2,3;-3,-1];
Jo = [3,3;-1,2];
F = J*inv(Jo);

% Part 4
w_tild = [1;0];
w = J*w_tild;
W = Jo*w_tild;

%Part 5
E_tild = .5*(J'*J-Jo'*Jo);

% Part 6 
e = inv(J)'*E_tild*inv(J);

%Part 7
E = inv(Jo)'*E_tild*inv(Jo);