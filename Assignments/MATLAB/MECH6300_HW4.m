%%% MECH 6300 - Homework 4

close all

% Problem 1

A = [0  1   0
     0  0   1
     -2 -4  -3];
B = [1  0
     0  1
     -1 1];
C = [0  1   -1
     1  2   0];
D = 0;

sys = ss(A,B,C,D)

% Matrix Inversions
syms s
sI_A = (s * eye(3) - A);
poles = factor(det(sI_A),s);
sI_A_inv = sI_A^-1


% Problem 2 and 3
A = [1  0   -1
    -16 -2  7
    0   -1  -2];
B = [1;0;1];
C = [1  -1  0];
D = 0;

sys = ss(A,B,C,D);
[csys,T] = canon(sys)

P_inv = [-1 -3  1
        2   20  0
        -2  -4  8];
P = inv(P_inv);

A_hat = P*A*P_inv;
B_hat = P*B;
C_hat = C*P_inv;
D_hat = D;

sys2 = ss(A_hat,B_hat,C_hat,D_hat)

figure()
step(sys)
hold on
step(csys)
step(sys2)
legend()
hold off


%Problem 4 Verification
A = [1  0   1   1
    0   1   0   0
    0   0   1   -1
    0   0   0   1];

s*eye(4)-A;
adjoint(s*eye(4)-A);
inv(s*eye(4)-A)
e_At = ilaplace(inv(s*eye(4)-A))
exp(A*t)


%Problem 6
A1 = [2 -1  2
      0 -1  -1
      0 0   1];
B1 = [1; 1; 0];
C1 = [1 -1 0];
D1 = 0;
sys1 = ss(A1,B1,C1,D1);
tf1 = tf(sys1)


A2 = [2 3   2
      0 -1  1
      0 0   -1];
B2 = [1; 1; 0];
C2 = [1 -1  0];
D2 = 0;
sys2 = ss(A2,B2,C2,D2);
tf2 = tf(sys2)


