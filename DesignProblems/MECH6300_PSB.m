%Problem 1
A = [   0   1   0       0;
        0   -8  -1.96   0;
        0   0   0       1;
        0   -8  11.76   0];
B = [0; 0.4; 0; -0.4];
C = [   1   0   0   0;
        0   0   1   0];
D = 0;
sys1 = ss(A,B,C,D)

%part 1b
syms t
exp(A*t)


%Problem 2
A = [   0   0   1       0;
        0   0   0       1;
        -40 -40 -400    0;
        40  -40 0       400];
B = [	0   0;
        0   0;
        2   0;
        0   2];
C = [   1   0   0   0;
        0   1   0   0];
D = 0;

sys2 = ss(A,B,C,D)

%part 2c
x0_1 = [0.5 1];
x0_2 = [-0.5 0.5];
Tf = 20;
Ts = 0.2;
t = linspace(0,Tf,Tf/Ts+1);
u_1a = 10 * gensig("square",1,Tf,Ts)-5;
u_1b = -3 * gensig("square",3,Tf,Ts) + 1.5;
u_1 = [u_1a, u_1b];
u_2a = 10 * gensig("square",2,Tf,Ts);
u_2b = 5 * gensig("square",5,Tf,Ts);
u_2 = [u_2a u_2b];

subplot(2,2,1)
lsim(sys2,u_1,t,x0_1)
subplot(2,2,2)
lsim(sys2,u_2,t,x0_1)
subplot(2,2,3)
lsim(sys2,u_1,t,x0_2)
subplot(2,2,4)
lsim(sys2,u_2,t,x0_2)

