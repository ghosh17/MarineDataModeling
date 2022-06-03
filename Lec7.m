%%clear, close, clc all function
clear all
close all
clc

load week3


%plot SST in x vs NO3

figure(1)
plot(SST, Y,'o');

%Chi Squared
a = [0.3418 0.2052];
Y1 = a(1) +a(2)*SST;

chiSq = sum((Y1-Y).^2/(std(Y).^2));

% %increase slope
% a(2) = 0.5;
% 
% Y1 = a(1) +a(2)*SST;


chiSq = sum((Y1-Y).^2/(std(Y).^2));
chiSwv = 1/(length(Y)-2)*chiSq;

hold on
plot(SST,Y1, '--')
hold off

%what to do if each term has diffrnt error
figure(2)
errorbar(SST,Y,Y_std,'o');%Y_std is Error estimate in the Y direction
hold on 
%some have large error bars some have small error bars
%plot(X0,Y1);%line we drew before error bars
%linfit can be given uncertainty in y
%need to download linfit

[a,sa,cov,r] = linfit(SST',Y',Y_std');

%hard way

ssy = 1./Y_std.^2;
S = sum(ssy);
Sx = sum(SST.*ssy);
Sxx = sum(SST.*SST.*ssy);
Sy = sum(Y.*ssy);
Sxy = sum(SST.*Y.*ssy);

%create A matrix
A=[S Sx;Sx Sxx];
rank(A);%Don't have to worry about rank deficeint  rank =2

b = [Sy;Sxy];
x = A\b;

%Can you linfit to the the same thing or can use regress if you dont have errors
%note it wants Y,X and assumes intercept of 0
[B,BINT,R,RINT] = regress(Y',[ones(size(SST))' SST']);




