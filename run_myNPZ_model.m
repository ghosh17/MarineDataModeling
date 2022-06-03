%%
close all
clear all
clc

%Define global variables
global p


p(1)= 0.4; %max phytoplankton growth rate (p(1)
p(2)  = .2; %Ksp p(2)
p(3) = .3;%Zooplankton grazing rate p(3)
p(4) = .2;%zooplankton p(4)
p(5) = .05; %m =Zooplankton mortality p(5)
p(6) = .2; %p(6)
p(7) = 2; %p(7)
p(8) = 20; %p(8)

%Or 

T = [0:0.1:400]';%how long to run. ODE solver needs column vector... 
n = length(T);
X0 = [0.8 .174 1]';%starting point for ODE solver. initial P,Z,N
[T,X] = ode15s('myNPZ_model',T,X0);
figure(1)
h = plot(T,X);

xlabel('time')
ylabel('[Nitrogen (mmol m-3]')
legend(h,'phyto','zoo','NO3')

%Need to define global variable so we can see what we defined outside the
%function inside the box