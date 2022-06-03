

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lotka Volterra

clear all
close all

figure(1)
N=[0:.1:4];
r01=.1;
K1=3;
plot(N, r01*N); hold on
plot(N, -r01/K1.*N.^2)
tmp=r01*N-r01/K1.*N.^2;
plot(N,  tmp, '--')
grid on
[I]=find(N==2);
H2=quiver(N(I),tmp(I),(N(I+3)-N(I)),(tmp(I+3)-tmp(I)),1,'k','linewidth',2,'maxheadsize',1);
[I]=find(N==4);
H2=quiver(N(I),tmp(I),(N(I-3)-N(I)),(tmp(I-3)-tmp(I)),1,'k','linewidth',2,'maxheadsize',1);
plot(K1, 0, 'o')

global a

%Case 1
close all
N0=[0.1 .1]';
r01=.1;  a(1)=r01;   % r01  N1 growth rate
K1=3;    a(2)=K1;    % K1   N1 carrying capacity
a12=2;  a(3)=a12;   % a12  effect of species 2 on species 1
r02=.1;  a(4)=r02;   % r02  N2 growth rate
K2=1;    a(5)=K2;    % K2   N2 carrying capacity
a21=1;  a(6)=a21;   % a21  effect of species 1 on species 2
[N]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21,N0);

% Case 2
close all
K1=1;    a(2)=K1;    % K1   N1 carrying capacity
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

% Case 3
close all
a21=3;  a(6)=a21;   % a21  effect of species 1 on species 2
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

N0=[1 1]';
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

N0=[.2 1]';
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

N0=[.1 .5]';
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

% Case 4
close all
a12=.1;  a(3)=a12;   % a12  effect of species 2 on species 1
K1=.2;    a(2)=K1;    % K1   N1 carrying capacity
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

N0=[.3 2]';
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

N0=[.02 1.4]';
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

r01=10;  a(1)=r01;   % r01  N1 growth rate
[N,t]=run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

(K1-a12*K2)/(K2-a21*K1)
N(end,1)/N(end,2)



%%% R vs K selected
close all
global a
N0=[.1 0.1]';
r01=.1;  a(1)=r01;   % r01  N1 growth rate
K1=5;    a(2)=K1;    % K1   N1 carrying capacity
a12=2;  a(3)=a12;   % a12  effect of species 2 on species 1
r02=.6;  a(4)=r02;   % r02  N2 growth rate
K2=2;    a(5)=K2;    % K2   N2 carrying capacity
a21=2;  a(6)=a21;   % a21  effect of species 1 on species 2
[N,t] = run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

r02=10;  a(4)=r02;   % r02  N2 growth rate
[N,t] = run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

r02=.1;  a(4)=r02;   % r02  N2 growth rate
[N,t] = run_LotkaVolterra (r01, K1, a12, r02, K2, a21, N0);

% How does R change the story?
% improve K1's ability to compete --> stable coexistance?

