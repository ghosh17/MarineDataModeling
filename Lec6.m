%%clear, close, clc all function
clear all
close all
clc


load week2

%plot(d4,e4,'*');
SPjk = sum((d4-mean(d4)).*(e4-mean(e4)));%Sum of products for j and k

SSj = sum((d4-mean(d4)).^2);

SSk = sum((e4-mean(e4)).^2);

rjk = SPjk/(sqrt(SSj*SSk)); %so weakly positively correlated

R = corrcoef(d4,e4); %3X3 if three variable, 2X2 if 2 variables

%NON Parametric test

%normplot(d4);
%normplot(e4);

%sort data for determing rank

[wq, rd4] = sort(d4);%first sorted vector, rd4 gives indexing 1st is smallest one then next smallest one.

%if wewanted to rearrange based on d4 sorted

sortedE4 = e4(rd4);%e4 is sorted based on smallest to largest of d4. The code makes new entry based on indexing retrived from sort d4. SO it is reordering


%instead of the sort(d4) we can use tiedrank, and if ties it will give an
%average rank

[id4, wq] = tiedrank(d4);

%if you wan to tconvert indexing from sort to come out of tiderank, u have
%to use ranking
id40(rd4) = 1:numel(d4);
1:numel(d4);

[ie4, wq] = tiedrank(e4);
%rs = sum((id4-mean(id4)).*(ie4-mean(ie4)))/(sqrt(sum((id4-mean(id4)).^2).*sqrt(sum((ie4-mean(ie4)).^2)));

%built in function

[RHO1, PVAL1] = corr (d4', e4', 'type', 'Spearman');

[RHO2, PVAL2] = corr (d4', e4', 'type', 'Kendall');








