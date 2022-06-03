function [pfinal] = ecopz_costfxn_opt_fn(cstr)

global tobs xobs wobs

p0(1)=.5;    % phytoplankton growth rate
p0(2)=2;     % phytoplankton carrying capacity
p0(3)=1.5;   % zooplankton grazing (Lokta-Voltera)
p0(4)=0.1;   % zooplankton mortality

preal(1)=1;     % phytoplankton growth rate
preal(2)=3;     % phytoplankton carrying capacity
preal(3)=2;     % zooplankton grazing (Lokta-Voltera)
preal(4)=0.3;   % zooplankton mortality

% look at results for first guess parameters
T=[0:0.1:100]';      % use a fixed time-scale to compare cases
n=length(T);
X0=[0.1 0.1]';

options = optimset('TolFun',10^-3); %let's not be too greedy
options = optimset('LargeScale','off'); % we don't have gradient information
%options = optimset('Display','iter'); %in case we want to see it all
[pfinal cost_final]=fminunc('ecopz_costfxn2',p0,options);


[preal' pfinal']

p=pfinal;
[T,X]=ode15s('ecopz',T,X0);
errorbar(tobs,xobs(:,1),wobs(:,1),'ko'); hold on
errorbar(tobs,xobs(:,2),wobs(:,2),'kv'); hold on
plot(T,X(:,1),'-','color', cstr);
plot(T,X(:,2),'-.','color', cstr);
ylabel('Plankton (mmol m^{-3})','fontsize',16);
xlabel('Time (days)','fontsize',16);
title('Final Optimized Fit','fontsize',16);
legend('Phytoplankton','Zooplankton','Phyto data','Zoop data','location','northeast');
legend boxoff;
set(gca,'fontsize',12);