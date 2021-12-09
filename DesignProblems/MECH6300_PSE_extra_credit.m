% MECH 6300 - Problem Set E
clear
close all

% Inverted Pendulum
m = 0.2;
M = 1;
l = 1;
g = 9.8;
k = 2;
R = 50;
p_obsv_ideal2 = 0.01;

A = [0  1   0   0;
    0   (-k^2)/(M*R*p_obsv_ideal2^2)    (-m*g)/(M)  0;
    0   0   0   1;
    0   (k^2)/(M*l*R*p_obsv_ideal2^2)   ((m+M)*g)/(M*l) 0];
B = [0;
    k/(M*R*p_obsv_ideal2);
    0;
    -k/(M*l*R*p_obsv_ideal2)];
C = [1  0   0   0];
D = 0;

sys = ss(A,B,C,D)

H = tf(sys)
H_zpk = zpk(H)

u = ctrb(A,B);
rank_u = rank(u)
v = obsv(A,C);
rank_v = rank(v)

sys_ctrbl = canon(sys,'companion')

% State-Feedback Design (from PSD)-----------------------------------
sys_feedback = ss(A,B,eye(4),0);
p_f_ideal = [-4,-2+j*2*sqrt(3),-2-j*2*sqrt(3),-800];

K = place(A,B,p_f_ideal)

newPoles_feedback = pole(feedback(sys_feedback,K))

% Simulation ---- Simulink Model
x0 = [5,0,0.2,0];

% Additional Speed Up
p_f_ideal_2 = [-40,-20+j*20*sqrt(3),-20-j*20*sqrt(3),-800];

K_2 = place(A,B,p_f_ideal_2);

newPoles_feedback2 = pole(feedback(sys_feedback,K_2));
% -----------------------------------------------------------------------

% Full-Order Observer
syms s
Delta_s_full = (s/5)^4 + 2.613 * (s/5)^3 + (2 + sqrt(2)) * (s/5)^2 + 2.613 * (s/5) + 1;
p_obsv_ideal = double(root(Delta_s_full,s));

figure()
scatter(real(p_obsv_ideal),imag(p_obsv_ideal))
xmin = -5; xmax = 0.5;
ymin = -5; ymax = -ymin;
xlim([xmin xmax])
ylim([ymin ymax])
line([0,0],ylim)
line(xlim,[0,0])
rectangle('Position',[-5,-5,10,10],'Curvature',[1 1])
title('Full-Order Observer Poles')
xlabel('Real')
ylabel('Imag')

close all

F = blkdiag([real(p_obsv_ideal(1)),imag(p_obsv_ideal(1)); ...
             -imag(p_obsv_ideal(1)),real(p_obsv_ideal(1))], ...
            [real(p_obsv_ideal(3)),imag(p_obsv_ideal(3)); ...
             -imag(p_obsv_ideal(3)),real(p_obsv_ideal(3))])
G = [1;1;1;1]
rank_ctrb_FG = rank(ctrb(F,G))

T = lyap(A,-F,-G*C)

det_T = det(T)

H = T * B

newPoles_obsv = eig(F)


% Reduced-Order Observer
syms s
Delta_s_reduced = (s/5)^3 + 2 * (s/5)^2 + 2 * (s/5) + 1;
p_obsv_ideal2 = double(root(Delta_s_reduced,s))

figure()
scatter(real(p_obsv_ideal2),imag(p_obsv_ideal2))
xmin = -6; xmax = 0.5;
ymin = -5; ymax = -ymin;
xlim([xmin xmax])
ylim([ymin ymax])
line([0,0],ylim)
line(xlim,[0,0])
rectangle('Position',[-5,-5,10,10],'Curvature',[1 1])
title('Reduced-Order Observer Poles')
xlabel('Real')
ylabel('Imag')

close all

C2 = [1,0,0,0];

F2 = blkdiag([real(p_obsv_ideal(1)),imag(p_obsv_ideal(1)); ...
             -imag(p_obsv_ideal(1)),real(p_obsv_ideal(1))], ...
            -5)
G2 = [1;1;1]
rank_ctrb_FG = rank(ctrb(F2,G2))

T2 = lyap(-F2,A,-G2*C2)

P2 = [C2;T2]

P2_inv = inv(P2)

H2 = T2 * B

newPoles_obsv = eig(F2)



close all


% Symbolic Version ------------------------------------------------------
syms s
s_I_A_inv = inv(s * eye(4) - A);
s_I_A_det = 25*det(s * eye(4) - A)

charPoly = factor(det(s * eye(4) - A),'FactorMode', 'real');

stm_no_den = simplify(25*(s_I_A_inv.* det(s*eye(4)-A)));

alphas = sym2poly(s_I_A_det/25)

u_bar_inv = [alphas(4),alphas(3),alphas(2),alphas(1);
            alphas(3),alphas(2),alphas(1),0;
            alphas(2),alphas(1),0,0;
            alphas(1),0,0,0]
P_inv = u * u_bar_inv
P = inv(P_inv)

A_bar = P * A * P_inv;
B_bar = P * B;
C_bar =     C * P_inv;
D_bar = D;

sys_conical = ss(A_bar,B_bar,C_bar,D_bar)

% State Feedback
alphas_bar_feedback = poly(p_f_ideal)

K_bar = fliplr(alphas(2:5)) - fliplr(alphas_bar_feedback(2:5))
K_sym = K_bar * P

newPoles_sym = pole(feedback(sys_feedback,-K_sym))



% Problem 3 (Extra Credit)---------------------------------------------
disp('Problem 3 (Extra Credit ----------------------------------')
P = sys;
L = -G;
K = -K;
syms s
G_sym = K * inv(s * eye(4) - (A + B*K + L*C)) * L


[symNum,symDen] = numden(G_sym(1)); %Get num and den of Symbolic TF
TFnum = sym2poly(symNum);    %Convert Symbolic num to polynomial
TFden = sym2poly(symDen);    %Convert Symbolic den to polynomial
G = tf(TFnum,TFden)


G_sys = ss(A+B*K+L*C,L,K,0)
G_zpk = zpk(G_sys)

figure()
rlocus(G*P)

figure()
rlocus(G*P)
xlim([-10,10])
ylim([-15,15])

T = 1 + G * P

T_zeros = zero(T)
T_poles = pole(T)


[Num,Den] = tfdata(G*P,'v');
syms s
sys_syms=poly2sym(Num,s)/poly2sym(Den,s);

syms a
T = 1+a*sys_syms;
roots = solve(T==0,a);

r = solve(roots==a,s)

a = 1.58
roots_at_a = double(subs(r))
