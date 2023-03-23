function [tout,yout] = EulerSolver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE ydot = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0=[a; b; c] is a column vector that contain the
% initial values y(0)=y0.
% OUTPUT:tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution
% at different grid points.

%TODO
%There is definietly some problem with yout, not too sure if it's correct.
%Run some more checks on it

tout = t(:);
N = length(t) - 1;
h = tout(2) - tout(1);
yout(1:3,1) = y0;
for n = 1:N
    tn = tout(n);
    yout(1:3,n+1) = yout(1:3,n)+h*f(tn,yout(1:3,n));
end