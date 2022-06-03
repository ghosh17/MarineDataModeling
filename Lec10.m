%%clear, close, clc all function
clear all
close all
clc

load week4

DATA = colstd(DATA0);

R = DATA' * DATA / (length(DATA) -1); %be careful use N or M don'r use leght(DATA) coz if M>N it's use M and vice versa

R1 = cov(DATA);

[V, Lambda] = eig(R); %lambda gives wright for eigen vactor.. Use this to find the axis that explains max var

%want to switch columns so that eigen vector with highest weight first and
%so on.

[lambda, ilambda] = sort(diag(Lambda), 'desc');%ilambda is index

%if you wnd up with very very small labmbda you might want to throw out.

V = V(:,ilambda);

Lambda = diag(lambda); 


S = sqrt(Lambda);

[U0,stmp,V0] = svd(DATA); %sign in PCA is meaningless coz you can flip around the axis. Can always switch it if it makes more sense

%percent of variance explained from PC axis comes from lambda
%trace sums along diagonal. trace(lambda) = sum(diag(lambda))
PoV = 100*diag(Lambda)/trace(Lambda); %Percent of variance. If PC1 axis is low that means you are lookng at noise.
sum(PoV);

%PoR = 100*diag(R)/trace(R); 

Ar = V*S;

figure(1)

h1 = plot(Ar(:,1), 'o-');
hold on
h2 = plot(Ar(:,2), '^-');
h3 = plot(Ar(:,3), 's-');
set(gca, 'xtick', [1 2 3]);
set(gca, 'xticklabel', [{'Var1'} {'Var2'} {'Var3'}]);
legend([h1 h2 h3], 'PC1', 'PC2', 'PC3');
%Here we see that variable 2 and 3 are highly weighted on PC1. 2 & 3 are
%positively related to PC1 and variable 1 is negatively related to PC1.


figure(2)
h1 = plot(V(:,1), 'o-');

hold on
h2 = plot(V(:,2), '^-');
h3 = plot(V(:,3), 's-');

set(gca, 'xtick', [1 2 3]);
set(gca, 'xticklabel', [{'Var1'} {'Var2'} {'Var3'}]);
legend([h1 h2 h3], 'PC1', 'PC2', 'PC3');
%Compared to figure 1 this doesnt tell us how much PC1 relative to PC2
%spread. 

Sr = DATA*V;%projects dataset onto new axis

figure(3)

h1 = plot(Sr(:,1), Sr(:,2), 'o');
hold on
xlabel('PC1');
ylabel('PC2');

%four samples so four points

%if your data has variability due to var 1 as like 1 but other varibale
%close to 0 then you can say var 1 driving variance. 

%Built in PCA
[COEFF,SCORE] = pca(DATA);%COEFF is V, SCORE is Sr

