rho = 28;
sigma = 10;
beta = 8/3;

f = @(t,x) lorenz(t,x,sigma,rho,beta);

k= 1;
h= 10^(-k);
tfinal=1;
t = [0:h:tfinal];

y0=[-1;3;4];

method = input('Available Methods:\n[1]:Euler\n[2]:Runge-Kutta order 4\n[3]:Implicit Runge-Kutta order 4\nChoose a method: ');
while ~((method ==1) || (method ==2) || (method ==3))
    method = input('Please enter 1, 2 or 3:');
end

if (method==1)
    [tout,Y] = EulerSolver(f,t,y0);
end
if (method==2)
    [tout,Y] = RK4Solver(f,t,y0);
end
if (method==3)
    [tout,Y] = IRK4Solver(f,t,y0);
end

options = odeset('RelTol',3.1e-14,'AbsTol',1e-16);
[toutm,Ym] = ode45(f,t,y0,options);

error = max(max(abs(Y-Ym')));

% 3d plot
plot3(Ym(:,1),Ym(:,2),Ym(:,3));
