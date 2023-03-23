function y = lorenz(t,x,sigma,rho,beta)
% INPUT: t is a a real value indicating time
%        x is a column vector of size 3 x 1
%        sigma, rho, beta are parameters of the Lorenz
%        equations
% OUPUT: y is a column vector of size 3 x 1 that gives
%        the right hand side of the Lorenz equations

if (length(x)~=3)
    error('x must contain 3 input values')
end


dxdt = sigma*(x(2)-x(1));
dydt = x(1)*(rho - x(3))-x(2);
dzdt = x(1)*x(2)-beta*x(3);
y = [dxdt ; dydt; dzdt];