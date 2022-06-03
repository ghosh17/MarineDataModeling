%%clear, close, clc all function
clear all
close all
clc

load week2

figure(1)
plot(x1,'r*')
hold on
plot(x2,'b.')

figure(2)
subplot(2,1,1) %2 rows, 1 cols. plot in the first plot
histogram(x1)
figure(2)
subplot(2,1,2)
histogram(x2)
figure(2)

[N1, E1, bin1] = histcounts(x1,10); %Number of data points from x1 that falls into each of the 10 bins, E1 is edges of bin, bin1 which bin each datapoint went into
size(N1);
bin1(1:5);
sum(N1); %100 as you're adding up all

[N2, E2, bin2] = histcounts(x2,10);


P1=N1/sum(N1);%what percentage of your data is in each bin

P2=N2/sum(N2);

%E1,E2 gives edges but the middle point is more helpful to plot
B1(1)=E1(1)+(E1(2)-E1(1))/2;
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;%Calculate center of the bins from the edges
B2=E2(1:end-1)+(E2(2:end)-E2(1:end-1))/2;

figure(3)
subplot(2,1,1)
bar(B1, P1)
hold on
subplot(2,1,2)
bar(B2,P2)
hold on

% you can make matlab calculate P1 (normalization probablility)
[P1, E1] = histcounts(x1,10, 'Normalization', 'probability');
subplot(2,1,1)
bar(B1,P1)
hold on
figure(3)

%%
close all
load week2
figure
plot(X1, '.')
size(X1); %million points!

%Sample 10 points randomly from million point dataset what does my dataset
%look like

[P1, E1] = histcounts(X1(1:10),bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;

figure(3)
subplot(2,1,1)
bar(B1,P1,'FaceColor', 'none', 'EdgeColor', 'b')
hold on
figure(3)

[P1, E1] = histcounts(X1(1:20),bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;
bar(B1,P1,'FaceColor', 'k')
figure(3)

%from going from 10 to 20 numbers we look more normal

[P1, E1] = histcounts(X1(1:30),bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;
bar(B1,P1,'FaceColor', 'm')
figure(3)

%20 looks less normal than 10

[P1, E1] = histcounts(X1(1:40),bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;
bar(B1,P1,'FaceColor', 'y')
figure(3)

[P1, E1] = histcounts(X1(1:50),bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;
bar(B1,P1,'FaceColor', 'g')
figure(3)

[P1, E1] = histcounts(X1(1:100),bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;
bar(B1,P1,'FaceColor', 'c')
figure(3)


%Whole million points
[P1, E1] = histcounts(X1,bin1, 'Normalization', 'probability');
B1=E1(1:end-1)+(E1(2:end)-E1(1:end-1))/2;
bar(B1,P1,'FaceColor', 'b')
figure(3)
%nice normal distribution. 
%As you sample more and more you get better at what the mean is ...

%X2 data

[P2, E2] = histcounts(X2,bin2, 'Normalization', 'probability');
B2=E2(1:end-1)+(E2(2:end)-E2(1:end-1))/2;
subplot(2,1,2)
bar(B2,P2,'FaceColor', 'b')
figure(3)
%Clearly not normally distributed


%Mean
mX1 = 1/length(X1) * sum(X1);
mX1 = mean(X1);
mX2 = mean(X2);

%Stdev
%sX1 = sqrt(1/length(X1)-1) * sum((X1-mean(X1)).^2);%wrote wrong
sX1 = std(X1);
sX2 = std(X2);


%varience
var(X1);

%Skewness
skewness(X1);

skewness(X2);

%

figure(11)
histogram(X1(1:200))
figure(11)
%what is the probablity that I measured a datapoint between 21.5 and 22 ie
%that bin. That is 25/200 fell into that bin so prob 0.1250. Equation 2.8
%probablitliy equation
X=14:28;
P0=1/((sX1)*sqrt(2*pi));
P = P0*exp(-(X-mX1).^2/(2*sX1.^2));
plot(X,P)


X=14:0.1:28;%taking this values
P0 = 1/((sX1)*sqrt(2*pi));
P = P0*exp(-(X-mX1).^2/(2*sX1.^2));
hold on
plot(X,P)
figure(1)
%matlab pdf: this is the matlab function fot the above plots
x1pdf = normpdf(X, mX1, sX1);
figure(1)
plot(X, x1pdf, '-.')
figure(1)
figure(3)
hold on
plot(X, x1pdf, 'k--')
figure(3)
subplot(2,1,1)


%%
%same thing continuation next class
close all
load week2

figure
plot(X1(1:200), '.')
figure(1)

%given dataset what P() that I made measurement of 21 - 21.5

figure(2)
ibin = 18:0.5:28;
histogram(X1, ibin)
figure(2)
%What % fell into bin
%43891/1000000;
%when we sub sampled it didn't look as gaussian.
%but these are so many points
%but we can use statistics


histogram(X1, ibin, 'Normalization', 'probability')

figure (2)

mX1 = mean(X1);

sX1 = std(X1);

x1pdf = normpdf(ibin(1):0.1:ibin(end), mX1, sX1);

figure(2)
hold on
plot(ibin(1):0.1:ibin(end), x1pdf)
figure(2)
%area under the curve = 1 
%By the mean point you have encountered half your data
%What fraction of your data have you encountered by 21

%Use cumulative dist: it adds aup probablility start at 0 

%hard way
x1cdf(1) = x1pdf(1);
for i=2:length(x1pdf)
    x1cdf(i) = x1cdf(i-1) + x1pdf(i)*0.1;%if u used a scale of 1 then dont need 0.1
end

plot(ibin(1):0.1:ibin(end), x1cdf)
figure(2)

%now you can read of cdf and see p <0.2 , p<0.3 etc

%easy way to do above

x1cpdf = normcdf(ibin(1):0.1:ibin(end), mX1, sX1);
plot(ibin(1):0.1:ibin(end), x1cpdf, '--')
figure(2)


%Check to see if normally distributed
figure(4)
normplot(X2)
%Dashed line is what you would expect for normal distribution ie linear
%distribution. If tou don't have normal data then it doen't fall on dashed
%line. This impacts which statistical test we can use


%What is the probablity of 21.5 given my mean and standard deviation
test1 = normcdf(21.5, mX1, sX1);
%What is the probablity of 21 given my mean and standard deviation
test2 = normcdf(21, mX1, sX1);
%What is the probability to fall between the two of them?
pBin = test1 - test2;
%This is what happens with t test

figure(6)
histogram(X2)
%This is not normally distributed. So we can't use operations for normal
%distribution 
figure(6)
subplot(2,1,1)
histogram(X2)
%using central limit therom. If you bin your data you will go from non
%normal to normal... read in book!!!
j = 1;%counter
%you are subsampling data. 

for i=1:10:length(X2) %want to bin in 10
    C2(j) = mean(X2(i:i+9));%average of those 10
    j = j + 1;
end
subplot(2,1,2)
histogram(C2)
figure(6)

normplot(C2);
%More normal

%So we went from non normal to normal distribution by binning.
    
    





























