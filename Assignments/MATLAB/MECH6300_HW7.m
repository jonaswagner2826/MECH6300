% MECH 6300 - HW 7

% Problem 1/2/3
n = 3;
A = [-1 -2  -3;
    0   -1  3;
    1   0   -1];
B = [3; 0; 2];
C = [3  3   0];

syms s

s_I_A_inv = inv(s * eye(n) - A);

charPoly = factor(det(s * eye(n) - A))%,'FactorMode', 'real')


% Problem 4
syms k1 k2
eq1 = 5 == -4 -3 * k1 -2 * k2;
eq2 = 6 == 11 + k1 + 2 * k2 + 6 * k1 * k2;

[k1,k2] = solve([eq1,eq2],[k1,k2]);
k1 = double(k1(1))
k2 = double(k2(1))


% Problem 5
A = blkdiag([2,1;0,2],-1,-1);
B = [0; 1; 1; 1];

% Part a
K = [1 1 1 1]
jordan(A + B*K)
K = [-1 3 5 1]
jordan(A + B*K)