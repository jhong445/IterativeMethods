function [tout,yout] = IRK4Solver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
%        the right-hand side of the ODE ydot = f(t,y)
%        t =[t0 t1 ... tfinal] is a vector of grid points
%           with length N
%        y0=[a; b; c] is a column vector that contain the
%        initial values y(0) = y0.
% OUTPUT:tout is a column vector of grid points.
%        yout is an 3 x N matrix containing the solution
%        at different grid points.


% The order 4 IRK table is 
% 1/2 - sqrt(3)/6 |     1/4           1/4 - sqrt(3)/6
% 1/2 + sqrt(3)/6 | 1/4 + sqrt(3)/6        1/4
%                 |----------------------------
%                 |   1/2                    1/2

% This leads to the following equations

% E_1 = y_n + h(1/4 * f(tn + (1/2 - sqrt(3)/6)*h,E_1) + (1/4 - sqrt(3)/6) * f(tn + (1/2 + sqrt(3)/6)*h,E_2))
% E_2 = y_n + h((1/4 + sqrt(3)/6)*f(tn + (1/2 - sqrt(3)/6)*h,E_1) + (1/4)*f(tn + (1/2 + sqrt(3)/6)*h,E_2))

% y_n+1 = y_n + h/2(f(tn + (1/2 - sqrt(3)/6)*h,E_1) + f(tn + (1/2 + sqrt(3)/6)*h,E_2))


h = t(2)-t(1);
N = length(t);
yout = zeros(3,N);
yout(:,1) = y0;
tout=t;

E_1 = @(yn,tn,E1,E2) yn + h*(1/4 * f(tn + (1/2 - sqrt(3)/6)*h,E1) + (1/4 - sqrt(3)/6) * f(tn + (1/2 + sqrt(3)/6)*h,E2)) -E1; 
E_2 = @(yn,tn,E1,E2) yn + h*((1/4 + sqrt(3)/6)*f(tn + (1/2 - sqrt(3)/6)*h,E1) + (1/4)*f(tn + (1/2 + sqrt(3)/6)*h,E2)) -E2;


for i = 1:N-1
X1 = @(S_E1,S_E2) E_1(yout(:,i),t(i),S_E1,S_E2);
X2 = @(S_E1,S_E2) E_2(yout(:,i),t(i),S_E1,S_E2);
X = @(S) [X1(S(1:3),S(4:6));X2(S(1:3),S(4:6))];
x0 = [0;0;0;0;0;0];

options = optimset('Display','off','TolFun',1e-15);

x = fsolve(X,x0,options);

E1 = x(1:3);
E2 = x(4:6);

yout(:,i+1) = yout(:,i) + (h/2)*(f(t(i) + (1/2 - sqrt(3)/6)*h,E1) + f(t(i) + (1/2 + sqrt(3)/6)*h,E2));
end