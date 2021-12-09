% MECH 6300 - HW 8

% Problem 1/2 ----------------------------------------------------------
n = 4;
A = [0  1   0   0;
    0   0   1   0;
    -2  1   3   1;
    1   2   0   0];
B = [0  0;
    0   0;
    2   1;
    0   1];

lambda_0 = eig(A)

syms s

s_I_A_inv = inv(s * eye(n) - A)

charPoly = factor(det(s * eye(n) - A),'FactorMode', 'real')

U = ctrb(A,B)
rank(U)

p = [-1+j*2,-1-j*2,-2+3*j,-2-3*j]

% Lyap Method
F = blkdiag([-1,-2;2,-1],[-2,-3;3,-2])
K_hat = [eye(2),eye(2)]

obsv_rank = rank(obsv(F,K_hat))

T = lyap(A,-F,B*K_hat)

det_T = det(T)

K_lyap = K_hat * inv(T)

eig_A_BK_lyap = eig(A+B*K_lyap)

% Place Method
K_place = place(A,-B,p)

eig_A_BK_place = eig(A + B * K_place)
%------------------------------------------------------------------------

% Problem 3/4 --------------------------------------------------------
A = [-3 -1  -2;
    0   -2  2;
    1   0   -2];
B = [2;0;1];
C = [-2,-1,0];

V = obsv(A,C)
rank_V = rank(V)

% Full Order
F = diag([-2,-3,-4])
G = ones(3,1)

u = ctrb(F,G)
rank_ctrl_FG = rank(u)

T = lyap(-F,A,-G*C)
det_T = det(T)

H = T * B

L = inv(T) * G

% Reduced Order
F = diag([-2,-4])
G = [1;1]

u = ctrb(F,G)
rank_u_FG = rank(u)

T = lyap(-F,A,-G*C)

H = T*B


