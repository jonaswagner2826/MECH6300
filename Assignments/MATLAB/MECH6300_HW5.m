%%% MECH 6300 - Homework 5

% Problem 3
charEq = [1 6   13  2   4   1]
poles3 = roots(charEq)


% Problem 4
A = [-3 1   0
    0   -2   0
    0   0   0];
B = [1
    2
    0];
C = [1  4   2];
D = 0;
sys = ss(A,B,C,D);
zpk4 = zpk(sys)


% Problem 5
A = [2  3   2
    3   1   0
    2   0   2];
eig_5a = eig(A)
syms s
delta_s = det(s*eye(3) - A)

A = [0  0   1
    0   0   0
    1   0   2];
Eig_5b = eig(A)