%Fisher model



for n=1:1000
time=100;
x=NaN*ones(time,1);
y=NaN*ones(time,1);
dx=NaN*ones(time,1);
dy=NaN*ones(time,1);

x(1)=.6; dx(1)=0;
y(1)=.6; dy(1)=0;
sigma=-.05;            %0.5
q(1)=sqrt(x(1).^2 +y(1).^2);

for t=1:time
    
    dx0=randn(1);
    dy0=randn(1);
    a=1/sqrt(dx0^2+dy0^2);
    dx(t+1)=sigma*a*dx0;
    dy(t+1)=sigma*a*dy0;
    x(t+1)=min(max(x(t)+dx(t+1),-100),100);
    y(t+1)=min(max(y(t)+dy(t+1),-100),100);
    
    q(t+1)=sqrt(x(t+1).^2 +y(t+1).^2);
    
    if q(t+1)>q(t)
        x(t+1)=x(t);
        y(t+1)=y(t);
        q(t+1)=q(t);
    end
    
      if n==0 && t<100
        figure(2)
        subplot(2,1,1)
        plot(x(1:t+1), y(1:t+1), 'k-'); hold on
        plot([x(t+1)], [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off

        subplot(2,1,2)
        plot([t t+1], 1-[q(t) q(t+1)], 'k')
        hold on
        pause(.1)
    end
    
end



if n<30
figure(1)
plot(x,y, 'k')
hold on
plot(x(end), y(end), 'mo', 'MarkerFaceColor', 'm');
pause(0.1)

figure(2)
plot(1:time+1, 1-q, 'k')
hold on
pause(0.1)
end

xend(n)=x(end);
yend(n)=y(end);
[I]=find(q<sigma);
if isempty(I)
    tend(n)=NaN;
else
    tend(n)=I(1);
end

end

figure(1)
plot(xend, yend, 'mo', 'MarkerFaceColor', 'm');

figure(3)
subplot(2,2,[1 2])
plot(xend, yend, 'mo', 'MarkerFaceColor', 'm');

subplot(2,2,3)
hist(xend)

subplot(2,2,4)
hist(yend)

figure
hist(tend)





SIGMA=[0.01 0.1 0.2 0.6 1 2 5];
MYcolor=['m' 'b' 'c' 'g' 'y' 'r' 'm'];

for s=1:length(SIGMA)

for n=1:1000
time=100;
x=NaN*ones(time,1);
y=NaN*ones(time,1);
dx=NaN*ones(time,1);
dy=NaN*ones(time,1);

x(1)=-10; dx(1)=0;
y(1)=-10; dy(1)=0;
sigma=SIGMA(s);          
q(1)=sqrt(x(1).^2 +y(1).^2);

for t=1:time
    
    dx0=randn(1);
    dy0=randn(1);
    a=1/sqrt(dx0^2+dy0^2);
    dx(t+1)=sigma*a*dx0;
    dy(t+1)=sigma*a*dy0;
    x(t+1)=min(max(x(t)+dx(t+1),-100),100);
    y(t+1)=min(max(y(t)+dy(t+1),-100),100);
    
    q(t+1)=sqrt(x(t+1).^2 +y(t+1).^2);
    
    if q(t+1)>q(t)
        x(t+1)=x(t);
        y(t+1)=y(t);
        q(t+1)=q(t);
    end
    
      if n==0 && t<100
        figure(2)
        subplot(2,1,1)
        plot(x(1:t+1), y(1:t+1), 'k-'); hold on
        plot([x(t+1)], [y(t+1)], 'ko', 'MarkerFaceColor', 'k')
        hold off

        subplot(2,1,2)
        plot([t t+1], [q(t) q(t+1)], 'k')
        hold on
        pause(.1)
    end
    
end

if n<30
figure(1)
plot(x,y, 'k')
hold on
plot(x(end), y(end), 'o', 'MarkerEdgeColor', MYcolor(s), 'MarkerFaceColor', MYcolor(s));
pause(0.1)

figure(2)
plot(1:time+1, q, MYcolor(s))
hold on
pause(0.1)
end

xend(s,n)=x(end);
yend(s,n)=y(end);
[I]=find(q<sigma);
if isempty(I)
    tend(s,n)=NaN;
else
    tend(s,n)=I(1);
end
qend(s,n)=q(end);
end

figure(1)
for ss=1:s
    if ss==length(SIGMA)
        plot(xend(ss,:), yend(ss,:), 'o', 'MarkerEdgeColor', MYcolor(ss)); hold on
    else
        plot(xend(ss,:), yend(ss,:), 'o', 'MarkerEdgeColor', MYcolor(ss), 'MarkerFaceColor', MYcolor(ss)); hold on
    end
end
end

figure
plot(SIGMA, mean(qend,2), 'o-')

