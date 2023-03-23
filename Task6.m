%% 3D plot of ERK with h=10^-3
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

[tout,Y] = RK4Solver(f,t,y0);

plot3(Y(1,:),Y(2,:),Y(3,:));


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


% ############
% Error tables
% ############



%% Error table for Euler

rho = 28;
sigma = 10;
beta = 8/3;

f = @(t,x) lorenz(t,x,sigma,rho,beta);


y0=[-1;3;4];

k_list = [2:1:6];
h_list = 10.^(-(k_list));
errors_EU = zeros(length(k_list),1);

options = odeset('RelTol',3.1e-14,'AbsTol',1e-16);

tfinal = 1;
fprintf("Errors table for Euler Method\n")
fprintf("h         |   Max Error\n")
fprintf("-----------------\n")

for i = [1:1:length(k_list)]
    h= h_list(i);
    t = [0:h:tfinal];

    [tout,Y] = EulerSolver(f,t,y0);
    [toutm,Ym] = ode45(f,t,y0,options);
    errors_EU(i) = max(max(abs(Y-Ym'))); 
    fprintf("%.1d   |     %e\n",h,errors_EU(i))
end

rho = zeros(length(k_list)-1,1);

fprintf("\nStep    |      Rho\n")
for i = [1:1:length(k_list)-1]
    
    rho(i) = log(errors_EU(i)/errors_EU(i+1))/log(h_list(i)/h_list(i+1));
    fprintf("E%d->E%d  |   %f\n",i,i+1,rho(i));
end



%% Error table for ERK

rho = 28;
sigma = 10;
beta = 8/3;

f = @(t,x) lorenz(t,x,sigma,rho,beta);


y0=[-1;3;4];

k_list = [1:1:4];
h_list = 10.^(-(k_list));
errors_ERK = zeros(length(k_list),1);

options = odeset('RelTol',3.1e-14,'AbsTol',1e-16);

tfinal = 1;
fprintf("Errors table for ERK4 Method\n")
fprintf("h         |   Max Error\n")
fprintf("-----------------\n") 

for i = [1:1:length(k_list)]
    h= h_list(i);
    t = [0:h:tfinal];

    [tout,Y] = RK4Solver(f,t,y0);
    [toutm,Ym] = ode45(f,t,y0,options);
    errors_ERK(i) = max(max(abs(Y-Ym'))); 
    fprintf("%.1d   |     %e\n",h,errors_ERK(i))
end

rho = zeros(length(k_list)-1,1);
fprintf("\nStep    |      Rho\n")
for i = [1:1:length(k_list)-1]
    rho(i) = log(errors_ERK(i)/errors_ERK(i+1))/log(h_list(i)/h_list(i+1));
    fprintf("E%d->E%d  |   %f\n",i,i+1,rho(i));
end

%% Error table for IRK

rho = 28;
sigma = 10;
beta = 8/3;

f = @(t,x) lorenz(t,x,sigma,rho,beta);


y0=[-1;3;4];

k_list = [1:1:4];
h_list = 10.^(-(k_list));
errors_IRK = zeros(length(k_list),1);

options = odeset('RelTol',3.1e-14,'AbsTol',1e-16);

tfinal = 1;
fprintf("Errors table for IRK4 Method\n")
fprintf("h         |   Max Error\n")
fprintf("-----------------\n")

for i = [1:1:length(k_list)]
    h= h_list(i);
    t = [0:h:tfinal];

    [tout,Y] = IRK4Solver(f,t,y0);
    [toutm,Ym] = ode45(f,t,y0,options);
    errors_IRK(i) = max(max(abs(Y-Ym'))); 
    fprintf("%.1d   |     %e\n",h,errors_IRK(i))
end

rho = zeros(length(k_list)-1,1);
fprintf("\nStep    |      Rho\n")
for i = [1:1:length(k_list)-1]
    
    rho(i) = log(errors_IRK(i)/errors_IRK(i+1))/log(h_list(i)/h_list(i+1));
    fprintf("E%d->E%d  |   %f\n",i,i+1,rho(i));
end


