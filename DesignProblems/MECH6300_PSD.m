% MECH 6300 - Problem Set D
% Inverted Pendulum
m = 0.2;
M = 1;
l = 1;
g = 9.8;
k = 2;
R = 50;
r = 0.1;

A = [0  1   0   0;
    0   (-k^2)/(M*R*r^2)    (-m*g)/(M)  0;
    0   0   0   1;
    0   (k^2)/(M*l*R*r^2)   ((m+M)*g)/(M*l) 0];
B = [0;
    k/(M*R*r);
    0;
    -k/(M*l*R*r)];
C = eye(4);%[1  0   0   0;
    %0   0   1   0];
D = 0;

sys = ss(A,B,C,D)

% Original Design
idealPoles = [-4,-2+j*2*sqrt(3),-2-j*2*sqrt(3),-800];

K = place(A,B,idealPoles);

newPoles = pole(feedback(sys,K));


% Simulation ---- Simulink Model
x0 = [5,0,0.2,0];

% Additional Speed Up
disp('Decreasing the distance to the jw axis to speed up response')
idealPoles_2 = [-40,-20+j*20*sqrt(3),-20-j*20*sqrt(3),-800]

K_2 = place(A,B,idealPoles_2)

newPoles2 = pole(feedback(sys,K_2))



% Problem 2 verification...
A = [0,0,1,0;
    0,0,0,1;
    -49,40,-400,0;
    40,-40,0,-400];
B = [0;0;2;0];
C=eye(4);
D=0;
sys =ss(A,B,C,D);

idealPoles3 = [-1-j,-1+j,-100-j*100,-100+j*100];
K3 = place(A,B,idealPoles3);
newPoles3 = pole(feedback(sys,K3));
