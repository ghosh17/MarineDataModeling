%%clear, close, clc all function
clear all
close all
clc

load week3


figure(1)
plot(SST,U,'o');
hold on

[a,sa,cov,r] = linfit(SST',U',0);

a0 =a;

[a,sa,cov,r] = linfit(U',SST',0); %a[1] order first intercept second a[2] is slope

%Depending on which axis you choose as independent you get differnt lines
%(differnt intercepts)


[m,b,r,sm,sb,xbar,ybar] = lsqfitma(SST,U);%minimizes a1 and a2.Type 2. m slope, b intecept


U0 = a0(1)+a0(2)+X0;
U0 = (1/a0(2))*(X0-a0(1));
plot(X0,U0);

UU = (1/a(2))*(X0-a(1));
plot(X0,UU,'--')


U2  =m*X0 + b;

plot(X0,U2,'-','linewidth',3)%suppose to show how it fits (goes through) both type on 1 lines


m2=sqrt(a0(2)/1/(a(2)));
%%
%%clear, close, clc all function
clear all
close all
clc


load phyto_growth.dat
% First time, phytoplankton abundance
figure(6)
plot(phyto_growth(:,1),phyto_growth(:,2), 'o');
hold on

%con trol axis
c = axis; %saves x and y axis value in variable xmin xmax ymin ymax

xlabel('Day','fontsize',14);
ylabel('Phytoplankton Growth', 'fontsize', 14);

%What to do if you want to fit a equation that we can't do design matrix for DO nlleasqr

%nlleasqr: what you give it is the function and it's going to brute force find minima
%a1 a2. First need to generate funciton

%create first function

time = [1:25];

pin = [0 1e4 0.5];

out = myFunc(time,pin);

plot(time, out)

axis(c);

%Non-linear least squares.. pin is guess

[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2] = nlleasqr(phyto_growth(:,1),phyto_growth(:,2),pin,'myFunc');

plot(phyto_growth(:,1),f)

figure(7)

Phyto = phyto_growth(:,2)/1e6;
plot(phyto_growth(:,1), Phyto, 'o');
hold on
c=axis;

pin = [0 0.01 0.5];

[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2] = nlleasqr(phyto_growth(:,1),Phyto,pin,'myFunc');


plot(phyto_growth(:,1),f,'--')

figure(8)

plot(phyto_growth(:,1), log(Phyto), 'o')
hold on

[a,sa,cov,r] = linfit(phyto_growth(:,1),log(Phyto),0);

P = a(1)+a(2).*phyto_growth(:,1);

plot(phyto_growth(:,1), P,'-');


