%%clear, close, clc all function
clear all
close all
clc


%ANOVA
load week2;
%J1 is same horizon

%Each row is differnt repeat

%Change from row to column. Left to right different sample. up down
%differnt repeats
DATA = [J1' J2' J3' J4' J5'];
[M,N] = size(DATA);

%Hard way
%SSa
SSa = sum((sum(DATA).^2)./M) - 1/(N*M)*(sum(sum(DATA))).^2;
SSt = sum(sum((DATA).^2)) - 1/(N*M) * sum(sum(DATA.^2));
SSw = SSt-SSa;
MSa = SSa / (N-1);
MSw = SSw / (N*(M-1));

F = MSa / MSw;

x = 0:0.1:10;
P_F = fcdf(F, (N-1), N*(M-1));

Prob_F=1-P_F;

y = fpdf(x, (N-1),N*(M-1));
plot(x,y);

%Easy way to do

[P, ANOVATAB,STATS] = anova1(DATA); %plots 5 samples and ANOVA table


%Two-way ANOVA
SSb = sum((sum(DATA,2).^2)./N) - 1/(N*M)*(sum(sum(DATA))).^2;
%mean(DATA,2)%WHat is the mean of all replicates of samples 1,2 ... the col 2 changes direction of mean

SSw = SSt -(SSa+SSb);
MSb = SSb/(M-1);
MSw = SSw/((N-1)*(M-1));
F1 = MSa/MSw;
F2=MSb/MSw;

P_F1 = fcdf(F1, (N-1), (N-1)*(M-1));
1-P_F1;
P_F2 = fcdf(F2, (M-1), (N-1)*(M-1));

%easy way
[P, ANOVATAB,STATS] = anova2(DATA);

