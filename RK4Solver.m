function [tout,yout] = RK4Solver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
%        the right-hand side of the ODE ydot = f(t,y)2
%        t =[t0 t1 ... tfinal] is a vector of grid points
%           with length N
%        y0=[a; b; c] is a column vector that contain the
%        initial values y(0) = y0.
% OUTPUT:tout is a column vector of grid points.
%        yout is an 3 x N matrix containing the solution
%        at different grid points

% The RK4 tableau from lec notes is
% \bf{c}
% 0   | 
% 1/2 | 1/2
% 1/2 |  0  1/2
% 1   |  0   0   1
% -------------------------
%\bf{b} 1/6 1/3 1/3  1/6

% So, we have
% E_1 = y_n
% E_2 = y_n + h*(1/2)*f(t_n,E_1)
% E_3 = y_n + h*(0)*... + h*(1/2)*f(t_n+c_2*h,E_2)
% E_4 = y_n + h*(0)*... + h*(0)*... + h*(1)*f(t_n+c_3*h,E_3)

% Then y_n+1 = y_n + h*[b_1*f(t_n+c_1*h,E_1) + b_2*f(t_n+c_2*h,E_2) + b_3*f(t_n+c_3*h,E_3) + b_4*f(t_n+c_4*h,E_4)]
% Then y_n+1 = y_n + h*[1/6*f(t_n+0*h,E_1) + 1/3*f(t_n+1/2*h,E_2) + 1/3*f(t_n+1/2*h,E_3) + 1/6*f(t_n+1*h,E_4)]


% Do some checks on input size, etc.
% TODO
if (length(y0)~=3)
    error('y0 must contain 3 input values')
end

tout = t(:);
h = tout(2)-tout(1);
N = length(t);

i=1;

yout(1:3,1) = y0;

while(i<N) 
    % Compute E_1 through E_4
    E_1 = yout(1:3,i);
    E_2 = yout(1:3,i) + h*(1/2) * f(t(i),E_1);
    E_3 = yout(1:3,i) + h*(1/2)*f(t(i)+(1/2)*h,E_2);
    E_4 = yout(1:3,i) + h*(1)*f(t(i)+(1/2)*h,E_3);

    % compute y_{n+1}
    yout(1:3,i+1) = yout(1:3,i) + h*((1/6) *f(t(i),E_1) + (1/3)*f(t(i)+(1/2)*h,E_2) + (1/3)*f(t(i)+(1/2)*h,E_3) + (1/6)*f(t(i)+h,E_4));

    i=i+1;
end 