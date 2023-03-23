%% 3D plot of IRK with h=10^-3
% (rho,sigma,beta) = (14,10,8/3)
% tfinal=100

rho = 14;
sigma = 10;
beta = 8/3;

y0=[-1;3;4];

f = @(t,x) lorenz(t,x,sigma,rho,beta);

k=3;
h= 10^(-k);
tfinal=100;
t = [0:h:tfinal];

[tout,Y] = IRK4Solver(f,t,y0);

plot3(Y(1,:),Y(2,:),Y(3,:));