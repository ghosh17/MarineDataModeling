%%clear, close, clc all function
clear all
close all
clc

load week4

%DATA0 matrix

mean(DATA0);

DATA = colstd(DATA0);%col stdardize

mean(DATA)%Means of 0

std(DATA);%STD of 1

%Calculate covarience explicitely 

[N,M] = size(DATA0);

for j=1:M
    for k =1:M
        sDATA(j,k) = (sum((DATA(:,j)-mean(DATA(:,j))).*(DATA(:,k)-mean(DATA(:,k)))))/(N-1);
        %diagonals become one, becasue we standardized it
    end
end

%slightly easier way to calc covarience
R = DATA'*DATA/(N-1); %R is another wat to calculate sDATA. Easier way

%Even easier way 
R1 = cov(DATA);

%when data standardized covarience and correlation coeff sae thing
R2 = corrcoef(DATA);

%PCA analysis

figure(2)
plot3(DATA(:,1),DATA(:,2),DATA(:,3),'*');
hold on
grid on
xlabel('variable 1');
ylabel('variable 2');
zlabel('variable 3');

[V,Lambda] = eig(R);%R covarience matrix

size(V);
Lambda;

