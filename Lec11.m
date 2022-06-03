%%clear, close, clc all function
clear all
close all
clc

%PCA
load week4.mat

%Davis0 is a length of boxes x and y
plot(Davis0(:,1), Davis0(:,2), 'o');
hold on

%****Peform PCA analysis on this data***

%Step 1: Standardize/Normalize data
Davis=colstd(Davis0);

%Step 2: Covarience matrix
R = cov(Davis);
figure(6)
plot(Davis(:,1), Davis(:,2),'bo');
hold on

%Step 3: Eigen
[V,Lambda] = eig(R);
%Lambda sorted wrong way
lambda = diag(Lambda);
[I,J] = sort(lambda, 'descend');
V = V(:,J);
%Remake the big Lambda
Lambda = diag(lambda(J));

%Step 4: Sr = Davis.V or Sf = Davis.Ar
Ar = V*sqrt(Lambda); %Ar is 2X2 dimension here
plot([-Ar(1,1) Ar(1,1)], [-Ar(2,1) Ar(2,1)], 'r');%PC1
plot([-Ar(1,2) Ar(1,2)], [-Ar(2,2) Ar(2,2)], 'r')%PC2

%What percent is described in each axis
PoV = 100*diag(Lambda)/trace(Lambda);% no of eigen value/vector = no of things measured

Sr = Davis*V;%projection of our dataset on vector
figure(7)
plot(Sr(:,1), Sr(:,2), 'o')
hold on
xlabel('PC1')
ylabel('PC2')

%Sf
figure(8)
Sf = Davis*Ar;
plot(Sf(:,1), Sf(:,2), 'o')

%Easy way
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(Davis);% in this built in function EXPLAINED gives PoV%
plot(SCORE(:,1), -SCORE(:,2), 'r*') %the +ve -ve can be flipped

%%
%Real dataset
%%clear, close, clc all function
clear all
close all
clc

%PCA
load week4.mat

DATA = BATS(:,2:8);%chose not to include depth measure for PCA
DATA = colstd(DATA);

R = cov(DATA);
[V,Lambda] = eig(R);
lambda = diag(Lambda);
[I,J] = sort(lambda, 'descend');
V = V(:,J);
%Remake the big Lambda
Lambda = diag(lambda(J));%Eigen values

%Step 4: Sr = Davis.V or Sf = Davis.Ar
Ar = V*sqrt(Lambda); %Ar is 2X2 dimension here
%plot([-Ar(1,1) Ar(1,1)], [-Ar(2,1) Ar(2,1)], 'r');%PC1
%plot([-Ar(1,2) Ar(1,2)], [-Ar(2,2) Ar(2,2)], 'r')%PC2

%What percent is described in each axis
PoV = 100*diag(Lambda)/trace(Lambda);

Sr = DATA*V;
figure(9)
plot(Sr(:,1), Sr(:,2), 'o')
hold on
xlabel('PC1')
ylabel('PC2')
%This says there is lot in this data that's unaccounted for. 
%7 dimensions plotted on 2 variab;es

plot([0 Ar(1,1)], [0 Ar(1,2)], '-', 'linewidth', 3)% temperature values
hold on
scatter(Sr(:,1), Sr(:,2), 50, DATA(:,1), 'filled');%50 is the size of the dot. So you can have to size based on a variable. The DATA term is color
%colored by temperature. We see Temp strongly related to PC1
plot([0 Ar(2,1)], [0 Ar(2,2)], 'k-', 'linewidth', 3)%Salinity
plot([0 Ar(3,1)], [0 Ar(3,2)], 'm-', 'linewidth', 3)%O2
plot([0 Ar(4,1)], [0 Ar(4,2)], 'c-', 'linewidth', 3)%
%^ if you get a plot like this you should break up the points shallow v/s
%deep samples

%Now we'll plot with depth. 
figure(10)
plot(Sr(:,1), BATS(:,1), 'o');%PC1 v/s depth plot. Therefore we think depth is driving PC1 in our data. So it's a good predictor of data

%Trace example: See lecture code
