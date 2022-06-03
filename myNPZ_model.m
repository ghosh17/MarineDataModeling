function [dXdt] = myNPZ_model(t, x) %ODE needs two input (1)->t (2)-> state varibales
dXdt = zeros(3,1);%make sure the vector in correct direction. Indexing still same. Alternatively you can say dXdt = dXdt' at the end...
%Define variables
global p


%State variables: N,P,Z

%ODE solvers don't need to know how many variables your solving
%We care about N, P , Z and we'll put into a vector
%where x(1) = P, x(2) = Z, x(3) = N 
%


%Equations

mu = p(1) * (x(3)/(x(3)+p(2)));

dXdt(1) = mu*x(1) - p(3)*x(1)*x(2);%dPdt calc

dXdt(2) =p(4) * p(3)*x(1)*x(2) - p(5)*x(2);%dZdt calc

dXdt(3) = p(6) * (p(7)-x(3))/p(8) - mu*x(1);%dNdt calc

