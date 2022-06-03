function xdot = ecopz(t,x)
% A function m-file for use with ODE45 to solve a set of coupled ordinary
% differential equations for a simple ecosystem model
%
% This is the simplified version, MLD is kept constant.
%
global MLD p
xdot=zeros(2,1);

% p(1) phytoplankton growth rate
% p(2) phytoplankton carrying capacity
% p(3) zooplankton grazing (Lokta-Voltera)
% p(4) zooplankton mortality

xdot(1) = p(1)*(1. -x(1)/p(2))*x(1) - p(3)*x(1)*x(2);
xdot(2) = p(3)*x(1)*x(2) - p(4)*x(2);
