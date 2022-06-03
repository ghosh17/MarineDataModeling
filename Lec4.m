%%clear, close, clc all function
%%
clear all
close all
clc

load week2

%S1 is notrmalized so it takes your data transforms mean to 0 and std to 1
mX1 = mean(X1);
sX1 = std(X1);
S1 = (X1-mX1)./sX1;
mean(S1);
std(S1);

xx = [1:length(x1)];

figure(1)
plot(xx,x1,'r*')

%say you went out and sampled and you got 25.3496 x =70 point. Is that a
%wierd data point? So we create z bin

%say you went out and sampled and you got x1 dataset. Is that a
%wierd data ? So we create z bin

zbins = [-5:0.001:5];
zpdf = normpdf(zbins);
figure(6)
plot(zbins, zpdf);%gaussian curve plotted
hold on

%What is the prob of getting value of 25.34 call it 25

%Z = (25-mX1); %mX1 is parent mean
%find mean and std of x1 dataset
mx1 = mean(x1);
sx1 = std(x1);

Z = (mx1 - mX1)/(sx1*sqrt(1/length(x1)));

%Now we have to see where the Z value 1.0508 falls

hold on 
c=axis; %gives axis of plot

figure(6)
plot([Z Z], [c(3) c(4)], 'r')

normcdf(-1.96); %given x value what is area under the curve 
norminv(0.025); %inverse of normcdf.. given area what is x value
 

critZ = norminv(0.025);
plot([critZ critZ], [c(3) c(4)], 'b')
plot([-critZ -critZ], [c(3) c(4)], 'b')
figure(6)
%if you are inside blue lines you are good.

%if you are within the lines you cannot reject null hypothesis!

figure(8)
subplot(2,1,1)
plot(d1,'r*')
hold on
plot(d2,'b*')
plot(d3,'c*')
plot(d4,'m*')

%Do these datasets come from the same million number dataset. Can we
%distinguish from parent

subplot(2,1,2)
h3 = histogram(d3, 'Normalization', 'probability');
E1 = h3.BinEdges;%Whateveredges you use for h3 histogram save them. We will use for the rest of the datasets.
hold on
h1 = histogram(d1, E1, 'Normalization', 'probability');
h2 = histogram(d2, E1, 'Normalization', 'probability');
h4 = histogram(d4, E1, 'Normalization', 'probability');
legend([h1 h2 h3 h4], ['D1', 'D2', 'D3', 'D4']);

norminv(0.025); %conf of 5% both tails 2.5% then inv
Z = (mean(d1)-23)/(std(d1)*sqrt(1/length(d1)));
%Z within +-1.96 so cannot reject null hypothesis 
Z = (mean(d2)-23)/(std(d2)*sqrt(1/length(d2)));
%Z -20 not within +-1.96 so reject null hypothesis 
Z = (mean(d3)-23)/(std(d3)*sqrt(1/length(d3)));
Z = (mean(d4)-23)/(std(d4)*sqrt(1/length(d4)));



H = ttest(d1,23);%H =0 that means null hypothesis cannot be rejected
H = ttest(d2,23);%H =1 that means null hypothesis can be rejected
H = ttest(d3,23);
H = ttest(d4,23);

 
%Two tail test
%Two samples
[H,P, CI, STATS] = ttest2(d1,d2); %P value is the type 1 alpha error








    