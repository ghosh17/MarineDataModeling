%http://www.math.cmu.edu/~shlomo/courses/BioSystems/Lectures.html

%Simple for loop
for n=1:1000 %run exact same walk a 1000 times
time=1000;%run for 1000 time steps
y=NaN*ones(time,1);%initialize y
y(1)=0;%start at 0
sigma=1;%size of the step we are taking
for t=1:time-1
%Cycle over for loop 1-1000
    y(t+1) = y(t)+sigma *randn(1); %location at next time step is current step plus some random movement. Becasue gaussian will mre likely to pick near 0
    %randn goes positve and negative....
    if n==1 && t<100 %plotting first 100 for clarity 
        figure(1)
        plot(1:t+1,y(1:t+1), 'k-'); hold on %plotting and little line segments
        plot(t+1, [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off
        pause(.1)
    end
end
figure(2)
plot(y)
hold on
yend(n)=y(end);
end

figure(3)
hist(yend)

%%%%%%%%%
%%
%%Withtout Gaussian
%http://www.math.cmu.edu/~shlomo/courses/BioSystems/Lectures.html

%Simple for loop
for n=1:1000 %run exact same walk a 1000 times
time=1000;%run for 1000 time steps
y=NaN*ones(time,1);%initialize y
y(1)=0;%start at 0
sigma=1;%size of the step we are taking
for t=1:time-1
%Cycle over for loop 1-1000
    tmp = rand([2 1]);
    if tmp(1) > 0.5 %50% of the time I want to do this and vice versa... if you want 80% probability uou do 0.8
        itmp = tmp(2);
    else
        itmp = -tmp(2);
    end
    %Want to go from -1 to 1... can't have rand becasue its only positive
    y(t+1) = y(t)+sigma *randn(1); %location at next time step is current step plus some random movement. Becasue gaussian will mre likely to pick near 0
    %randn goes positve and negative....
    if n==1 && t<100 %plotting first 100 for clarity 
        figure(1)
        plot(1:t+1,y(1:t+1), 'k-'); hold on %plotting and little line segments
        plot(t+1, [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off
        pause(.1)
    end
end
figure(2)
plot(y)
hold on
yend(n)=y(end);
end

figure(3)
histogram(yend)


%%
clear yend xend x y

for n=1:10
time=1000;
x=NaN*ones(time,1);
y=NaN*ones(time,1);
z=NaN*ones(time,1);
x(1)=0;
y(1)=0;
z(1)=0;
sigma=1;
for t=1:time
%
    tmp = rand([2,1]);%determine sign for x and y
    [i] = find(tmp<0.5);
    tmp2 = ones(2,1);
    tmp2(i) = -1; %gets set to -1 if <0.5
    
    dx = sigma*tmp2(1)*rand(1);
    dy = sqrt(sigma^2 - dx^2)*tmp2(2);
    
    %new step
    x(t+1) = x(t) + dx;
    y(t+1) = y(t) + dy;
    %check step
    z(t+1) = sqrt((x(t+1)-x(t))^2 + (y(t+1)-y(t))^2);%All of them should be one becasue we didn't move in the z direction. Good way to check
    
    
    if n==1 && t<100
        figure(1)
        plot(x(1:t+1), y(1:t+1), 'k-'); hold on
        plot([x(t+1)], [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off
        pause(.1)
    end
end
plot(x,y)
pause(.1)
hold on
xend(n)=x(end);
yend(n)=y(end);
end

figure
hist(yend)

%%
    

%%% with sensing
%Bunch of bactiria introduced at 0,0 
%if they move forward and like it they will go forward, and if they don't like their step they will roll dice
%again.
clear yend xend x y

[X, Y] = meshgrid(0:100);
Q=.01*X+.01*Y;
figure (1)
h=pcolor(X,Y,Q);
set(h,'edgecolor','none')
hold on

%SAME LOOP FROM ABOVE
for n=1:100
time=1000;
x=NaN*ones(time,1);
y=NaN*ones(time,1);
dx=NaN*ones(time,1);
dy=NaN*ones(time,1);

x(1)=1; dx(1)=0;
y(1)=1; dy(1)=0;
sigma=1;
for t=1:time
    
    tmp=rand([2,1]);  % determine sign for x and y
    [I]=find(tmp<0.5);
    tmp2=ones(2,1);
    tmp2(I)=-1;
    dx=tmp2(1)*rand(1);
    dy=sqrt(sigma^2 - dx^2)*tmp2(2);
    x(t+1)= min(max(x(t) + dx,1),100);
    y(t+1)= min(max(y(t) + dy,1),100);
    q0=Q(floor(x(t)), floor(y(t)));
    q1=Q(floor(x(t+1)), floor(y(t+1)));
    if q0>q1
        x(t+1)=x(t);
        y(t+1)=y(t);
    end
    q(t+1)=Q(floor(x(t+1)), floor(y(t+1)));
     

      if n==1 && t<100
        figure(2)
        subplot(2,1,1)
        h=pcolor(X,Y,Q);
        set(h,'edgecolor','none')
        hold on
        plot(x(1:t+1), y(1:t+1), 'k-')
        plot([x(t+1)], [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off

        subplot(2,1,2)
        plot([t t+1], [q0 q1], 'k')
        hold on
        pause(.1)
    end
    
end
figure(1)
plot(x,y, 'k')
hold on
plot(x(end), y(end), 'mo', 'MarkerFaceColor', 'm')

figure(2)
plot(1:time+1, q, 'k')
hold on
[I]=find(q==max(q));
xend(n)=x(end);
yend(n)=y(end);
end
%%


%%% with chemotaxis

%Instead of sensing and rejecting, if the result is good I'll run for
%longer!

%SAME LOOP FROM ABOVE
clear yend xend x y

[X, Y] = meshgrid(0:100);
Q=.01*X+.01*Y;
figure (1)
h=pcolor(X,Y,Q);
set(h,'edgecolor','none')
hold on

%SAME LOOP FROM ABOVE
for n=1:100
time=1000;
x=NaN*ones(time,1);
y=NaN*ones(time,1);
dx=NaN*ones(time,1);
dy=NaN*ones(time,1);

x(1)=1; dx(1)=0;
y(1)=1; dy(1)=0;
sigma=1;
sigmaGO=10;
for t=1:time
    tmp=rand([2,1]);  % determine sign for x and y
    [I]=find(tmp<0.5);
    tmp2=ones(2,1);
    tmp2(I)=-1;
    dx(t+1)=tmp2(1)*rand(1);
    dy(t+1)=sqrt(sigma^2 - dx(t+1)^2)*tmp2(2);
    x(t+1)= min(max(x(t) + dx(t+1),1),100);
    y(t+1)= min(max(y(t) + dy(t+1),1),100);
    q0=Q(floor(x(t)), floor(y(t)));
    q1=Q(floor(x(t+1)), floor(y(t+1)));
    if q1>q0
        dx(t+1)=dx(t+1)*sigmaGO/sigma;
        dy(t+1)=dy(t+1)*sigmaGO/sigma;
        x(t+1)= min(max(x(t) + dx(t+1),1),100);
        y(t+1)= min(max(y(t) + dy(t+1),1),100);
    end
    q(t+1)=Q(floor(x(t+1)), floor(y(t+1)));

    
    
    if n==1 && t<100
        figure(2)
        subplot(2,1,1)
        h=pcolor(X,Y,Q);
        set(h,'edgecolor','none')
        hold on
        plot(x(1:t+1), y(1:t+1), 'k-')
        plot([x(t+1)], [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off

        subplot(2,1,2)
        plot([t t+1], [q0 q1], 'k')
        hold on
        pause(.1)
    end
    
end
figure(1)
plot(x,y, 'k')
hold on
plot(x(end), y(end), 'mo', 'MarkerFaceColor', 'm')

figure(3)
plot(1:time+1, q, 'k')
hold on
[I]=find(q==max(q));
t_end(n)=I(1);
end


%%

%%% with chemotaxis
%trying to hit the middle of the plot, more complicated gradient 
%THIS IS EXAMPLE NEEDED FOR PROBLEM SET

clear yend xend x y

[X, Y] = meshgrid(0:100);
Q=-cos(X/100*2*pi)-cos(Y/100*2*pi);
Q(Q<0)=0;
figure (1)
h=pcolor(X,Y,Q);
set(h,'edgecolor','none')
hold on

for n=1:100
time=1000;
x=NaN*ones(time,1);
y=NaN*ones(time,1);
dx=NaN*ones(time,1);
dy=NaN*ones(time,1);

x(1)=1; dx(1)=0;
y(1)=1; dy(1)=0;
sigma=1;    %5
sigmaGO=10;
for t=1:time
    tmp=rand([2,1]);  % determine sign for x and y
    [I]=find(tmp<0.5);
    tmp2=ones(2,1);
    tmp2(I)=-1;
    dx(t+1)=tmp2(1)*rand(1);
    dy(t+1)=sqrt(sigma^2 - dx(t+1)^2)*tmp2(2);
    x(t+1)= min(max(x(t) + dx(t+1),1),100);
    y(t+1)= min(max(y(t) + dy(t+1),1),100);
    q0=Q(floor(x(t)), floor(y(t)));
    q1=Q(floor(x(t+1)), floor(y(t+1)));
    if q1>q0
        dx(t+1)=dx(t+1)*sigmaGO/sigma;
        dy(t+1)=dy(t+1)*sigmaGO/sigma;
        x(t+1)= min(max(x(t) + dx(t+1),1),100);
        y(t+1)= min(max(y(t) + dy(t+1),1),100);
    end
    q(t+1)=Q(floor(x(t+1)), floor(y(t+1)));

    if n==1 %&& t<200
        figure(2)
        subplot(2,1,1)
        h=pcolor(X,Y,Q);
        set(h,'edgecolor','none')
        hold on
        plot(x(1:t+1), y(1:t+1), 'k-')
        plot([x(t+1)], [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off

        subplot(2,1,2)
        plot([t t+1], [q0 q1], 'k')
        hold on
        pause(.1)
    end
    
end
figure(1)
plot(x,y, 'k')
hold on
plot(x(end), y(end), 'mo', 'MarkerFaceColor', 'm')
pause(.1)

figure(3)
plot(1:time+1, q, 'k')
hold on
[I]=find(q==max(q));
if max(q)<1
    t_end(n)=NaN;
else
    t_end(n)=I(1);
end

end

[I]=find(isnan(t_end));
hist(t_end, [min(t_end):100:max(t_end)])
