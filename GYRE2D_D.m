function [c, iso, km, Oxysat, Delfin, LossTerm, Inventory, u, v, A, km2] = GYRE2D_D (J,a, km2,vect)


%GYRE2D_D
%isopycnal slab model with double gyres with an intensified western
%boundary
%increased diffusive mixing between two gyres at western boundary
%Southern boundary is restored to constant value
%Eastern, western and northern boundaries are treated as walls
%the consumption term J) and the fractionation term (a) need to be set
%when concentration drops below 2.5 (1% of sat) J=0

warning off 				% to take care of v4 vs. v5 warnings
tyr=3600*24*365.25; nx=50; ny=50;	% various constants
%dx=50000; dy=dx; 			% node spacing
%dx=25000; dy=dx; 			% node spacing
dx=60000; dy=dx; 			% node spacing


%Coordinates
tmp=-1:1/(nx/2):1; tmp2=1-.1*tmp.^2;

y=meshgrid([0:ny]/(ny));  x=meshgrid([0:nx]/(nx)); y=y';
y2=meshgrid([0:ny/2]/(ny/2), [0:nx]/(nx));  x2=meshgrid([0:nx]/(nx), [0:ny/2]/(ny/2)); y2=y2';

%Stommel gyre

%A=.0035;   %mean abs(u), mean abs(v) =.002
%A=.0021;    %mean vector (u,v) =.002
%A=0.0029;      %mean vector (u,v) =.002 for half field u(:,25:end), v(:,25:end) -western boundary excluded
%A=0.0143;      %mean vector (u,v) =.01 for half field u(:,25:end), v(:,25:end) -western boundary excluded
%A=0.0071;       %mean vector (u,v) =.005 for half field u(:,25:end), v(:,25:end) -western boundary excluded

A=1;
epsilon=.1;
lamda1=(-1+sqrt(1+4*(pi^2)*(epsilon^2)))/(2*epsilon);
lamda2=(-1-sqrt(1+4*(pi^2)*(epsilon^2)))/(2*epsilon);
c1=(1-exp(lamda2))/(exp(lamda2)-exp(lamda1));
c2=-(1+c1);


for j=1:nx+1
    for i=1:ny+1
        SF1(i,j)=A.*sin(pi*y(i,j)).*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
    end
%    for i=1:ny+1
%        SF2(i,j)=A.*sin(pi*y2(i,j)).*(c1.*exp(lamda1.*x2(i,j)) + c2.*exp(lamda2.*x2(i,j)) +1);
%    end
end

SF=[SF1;SF1(2:ny+1,:)];


u1=zeros(ny+1, nx+1); v1=u1;
u2=zeros(ny/2+1, nx+1); v2=u2;

for j=1:nx+1
    for i=1:ny+1
        yy = y(i,j)+1/(2*ny);
        xx = x(i,j)+1/(2*nx);
        u1(i,j)=pi.*cos(pi.*yy).*A.*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
        v1(i,j)=-sin(pi.*y(i,j)).*A.*(lamda1.*c1.*exp(lamda1.*xx)+ lamda2.*c2.*exp(lamda2.*xx));
        
        u2(i,j)=-pi.*cos(pi.*yy).*A.*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
        v2(i,j)=sin(pi.*y(i,j)).*A.*(lamda1.*c1.*exp(lamda1.*xx)+ lamda2.*c2.*exp(lamda2.*xx));
    end
    
%    for i=1:ny/2
%        yy = y2(i,j)+1/(2*(ny/2));
%        xx = x2(i,j)+1/(2*nx);
%        u2(i,j)=-pi.*cos(pi.*yy).*A.*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
%        v2(i,j)=sin(pi.*y(i,j)).*A.*(lamda1.*c1.*exp(lamda1.*xx)+ lamda2.*c2.*exp(lamda2.*xx));
%    end

end



u=[u1; u2(2:ny+1,:)];  v=[v1; v2(2:ny+1,:)];
tmp=sqrt(u.^2+v.^2);

A=vect/mean(mean(tmp(:,25:end)));

for j=1:nx+1
    for i=1:ny+1
        SF1(i,j)=A.*sin(pi*y(i,j)).*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
    end
%    for i=1:ny+1
%        SF2(i,j)=A.*sin(pi*y2(i,j)).*(c1.*exp(lamda1.*x2(i,j)) + c2.*exp(lamda2.*x2(i,j)) +1);
%    end
end

SF=[SF1;SF1(2:ny+1,:)];
figure(1)
subplot(2,2,1)
contour(SF)
title('Stream Function')

u1=zeros(ny+1, nx+1); v1=u1;
u2=zeros(ny/2+1, nx+1); v2=u2;

for j=1:nx+1
    for i=1:ny+1
        yy = y(i,j)+1/(2*ny);
        xx = x(i,j)+1/(2*nx);
        u1(i,j)=pi.*cos(pi.*yy).*A.*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
        v1(i,j)=-sin(pi.*y(i,j)).*A.*(lamda1.*c1.*exp(lamda1.*xx)+ lamda2.*c2.*exp(lamda2.*xx));
        
        u2(i,j)=-pi.*cos(pi.*yy).*A.*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
        v2(i,j)=sin(pi.*y(i,j)).*A.*(lamda1.*c1.*exp(lamda1.*xx)+ lamda2.*c2.*exp(lamda2.*xx));
    end
    
%    for i=1:ny/2
%        yy = y2(i,j)+1/(2*(ny/2));
%        xx = x2(i,j)+1/(2*nx);
%        u2(i,j)=-pi.*cos(pi.*yy).*A.*(c1.*exp(lamda1.*x(i,j)) + c2.*exp(lamda2.*x(i,j)) +1);
%        v2(i,j)=sin(pi.*y(i,j)).*A.*(lamda1.*c1.*exp(lamda1.*xx)+ lamda2.*c2.*exp(lamda2.*xx));
%    end

end



u=[u1; u2(2:ny+1,:)];  v=[v1; v2(2:ny+1,:)];
xx=[x; x(2:ny+1,:)];  yy=[y; y(2:ny+1,:)];

tmp=sqrt(u.^2+v.^2);

figure(1)
subplot(2,2,2)
contour(u)
title('u')
subplot(2,2,3)
contour(v)
title('v')

subplot(2,2,4)
ix=1:6:length(xx(:,1));
iy=1:6:length(yy(1,:));
quiver(u(ix,iy)*1000,v(ix,iy)*1000)
title('velocity')



%u=u1; v=v1;
%ny=ny+ny/2+1;
ny=ny+ny;

X=[0:nx-1]*dx/1000;Y=[0:ny-1]*dy/1000;	% distance vectors for plotting (in km)

km=1000;				% lateral diffusivity of 1000m2/s
%km2=7000;
vmax=max(max(v)); umax=max(max(u));		% maximum possible velocities

%dt=0.2*(dx/(umax+vmax+4*km2/dx));		% stable time step
dt=.5*(dx/(umax+vmax+4*km2/dx));		% stable time step

nt=20*tyr/dt;					% number of steps in 20 years
kx=km*ones(ny+1,nx+1); ky=kx;			% kx, ky over all domain

kx(40:70,1:25)=km2;
ky(40:70,1:25)=km2;


%OxyEquil=250*61.795850754258083/100;
OxyEquil=250;
delEquil=.75;
%delEquil=5.312;
startdel= ((delEquil/1000+1)*0.002)*OxyEquil;
c=OxyEquil*ones(ny,nx); cold=c;	iso=startdel*ones(ny,nx);	isoold=iso;		% initialize concentrations to zero


%concentration gradient
%c=meshgrid([1:2:ny*2],[1:nx]); c=c';
%cold=c;

%dye test
%c=zeros(ny,nx); 
%c(10,10)=1;
%cold=c;

		% now calculate weight matrices

wxm=zeros(ny,nx); wym=wxm; w0=wxm;

w0=1-kx(1:ny,1:nx)*dt/dx/dx -kx(1:ny,2:nx+1)*dt/dx/dx -ky(1:ny,1:nx)*dt/dy/dy- ky(2:ny+1,1:nx)*dt/dy/dy;
wxm=kx(1:ny,1:nx)*dt/dx/dx;
wxp=kx(1:ny,2:nx+1)*dt/dx/dx;
wym=ky(1:ny,1:nx)*dt/dy/dy;
wyp=ky(2:ny+1,1:nx)*dt/dy/dy;


for i=1:ny
    for j=1:nx
        if u(i,j)>=0 & u(i,j+1)>=0
            w0(i,j)=w0(i,j)-u(i,j+1)*dt/dx;
            wxm(i,j)=wxm(i,j)+u(i,j)*dt/dx;
        elseif u(i,j)<=0 & u(i,j+1)<=0
            w0(i,j)=w0(i,j)+u(i,j)*dt/dx;
            wxp(i,j)=wxp(i,j)-u(i,j+1)*dt/dx;
        elseif u(i,j)>0 & u(i,j+1)<0
            wxm(i,j)=wxm(i,j)+u(i,j)*dt/dx;
            wxp(i,j)=wxp(i,j)-u(i,j+1)*dt/dx;
        elseif u(i,j)<0 & u(i,j+1)>0
            w0(i,j)=w0(i,j)-u(i,j+1)*dt/dx+u(i,j)*dt/dx;
        end
        
        if v(i,j)>=0 & v(i+1,j)>=0
            w0(i,j)=w0(i,j)-v(i+1,j)*dt/dy;
            wym(i,j)=wym(i,j)+v(i,j)*dt/dy;
        elseif v(i,j)<=0 & v(i+1,j)<=0
            w0(i,j)=w0(i,j)+v(i,j)*dt/dy;
            wyp(i,j)=wyp(i,j)-v(i+1,j)*dt/dy;
        elseif v(i,j)>0 & v(i+1,j)<0
            wym(i,j)=wym(i,j)+v(i,j)*dt/dy;
            wyp(i,j)=wyp(i,j)-v(i+1,j)*dt/dy;
        elseif v(i,j)<0 & v(i+1,j)>0
            w0(i,j)=w0(i,j) - v(i+1,j)*dt/dy + v(i,j)*dt/dy;
        end
    end
end



%Boundary conditions
wxm(:,1)=0;
wxp(:,nx)=0;
wym(1,:)=0;
wyp(ny,:)=0;

w0(:,1)=w0(:,1)+ kx(1:ny,1)*dt/dx/dx;  %Western
w0(:,nx)=w0(:,nx)+ kx(1:ny,nx)*dt/dx/dx;  %Eastern
w0(1,:)=w0(1,:)+ ky(1,1:nx)*dt/dy/dy;  %South
w0(ny,:)=w0(ny,:)+ ky(ny,1:nx)*dt/dy/dy;  %North



        

J=J*1.026/tyr;  %convert umol/kg/yr to mmol/m3/sec
Jterm=J*dt*ones(ny,nx);
Jiso=J*dt*a*ones(ny,nx);


oldyear=1;				% our trusty year counter
%maxError= 0.002;
maxError= 0.02;

t=1;
r=1;
%%
    consum=Jterm;
    Isoconsum=Jiso;
    cold(1:2,:)=OxyEquil*ones(2,nx);		% set oxygen boundary conditions to 250
    isoold(1:2,:)=startdel*ones(2,nx);		% set O18 boundary condition to .75

    newyear=ceil(t*dt/tyr);		% time counter

    c(2:ny-1,2:nx-1) = w0(2:ny-1,2:nx-1).*cold(2:ny-1,2:nx-1) +...
	wxm(2:ny-1,2:nx-1).*cold(2:ny-1,1:nx-2) +...
	wxp(2:ny-1,2:nx-1).*cold(2:ny-1,3:nx) +...
	wym(2:ny-1,2:nx-1).*cold(1:ny-2,2:nx-1) +...
	wyp(2:ny-1,2:nx-1).*cold(3:ny,2:nx-1)+ consum(2:ny-1,2:nx-1);		% the main event


    iso(2:ny-1,2:nx-1) = w0(2:ny-1,2:nx-1).*isoold(2:ny-1,2:nx-1) +...
	wxm(2:ny-1,2:nx-1).*isoold(2:ny-1,1:nx-2) +...
	wxp(2:ny-1,2:nx-1).*isoold(2:ny-1,3:nx) +...
	wym(2:ny-1,2:nx-1).*isoold(1:ny-2,2:nx-1) +...
	wyp(2:ny-1,2:nx-1).*isoold(3:ny,2:nx-1)+ Isoconsum(2:ny-1,2:nx-1).*isoold(2:ny-1,2:nx-1)./cold(2:ny-1,2:nx-1);		% the main event

%Western Boundary
    c(2:ny-1,1) = w0(2:ny-1,1).*cold(2:ny-1,1) +...
	wxp(2:ny-1,1).*cold(2:ny-1,2) +...
	wym(2:ny-1,1).*cold(1:ny-2,1) +...
	wyp(2:ny-1,1).*cold(3:ny,  1)+ consum(2:ny-1,1);		% the main event

    iso(2:ny-1,1) = w0(2:ny-1,1).*isoold(2:ny-1,1) +...
	wxp(2:ny-1,1).*isoold(2:ny-1,2) +...
	wym(2:ny-1,1).*isoold(1:ny-2,1) +...
	wyp(2:ny-1,1).*isoold(3:ny,1)+ Isoconsum(2:ny-1,1).*isoold(2:ny-1,1)./cold(2:ny-1,1);		% the main event

%Eastern Boundary
    c(2:ny-1,nx) = w0(2:ny-1,nx).*cold(2:ny-1,nx) +...
	wxm(2:ny-1,nx).*cold(2:ny-1,nx-1) +...
	wym(2:ny-1,nx).*cold(1:ny-2,nx) +...
	wyp(2:ny-1,nx).*cold(3:ny,  nx)+ consum(2:ny-1,nx);		% the main event

    iso(2:ny-1,nx) = w0(2:ny-1,nx).*isoold(2:ny-1,nx) +...
	wxm(2:ny-1,nx).*isoold(2:ny-1,nx-1) +...
	wym(2:ny-1,nx).*isoold(1:ny-2,nx) +...
	wyp(2:ny-1,nx).*isoold(3:ny,nx)+ Isoconsum(2:ny-1,nx).*isoold(2:ny-1,nx)./cold(2:ny-1,nx);		% the main event

%Southern Boundary
    c(1,2:nx-1) = w0(1,2:nx-1).*cold(1,2:nx-1) +...
	wxm(1,2:nx-1).*cold(1,1:nx-2) +...
	wxp(1,2:nx-1).*cold(1,3:nx) +...
	wyp(1,2:nx-1).*cold(2,2:nx-1)+ consum(1,2:nx-1);		% the main event

    iso(1,2:nx-1) = w0(1,2:nx-1).*isoold(1,2:nx-1) +...
	wxm(1,2:nx-1).*isoold(1,1:nx-2) +...
	wxp(1,2:nx-1).*isoold(1,3:nx) +...
	wyp(1,2:nx-1).*isoold(2,2:nx-1)+ Isoconsum(1,2:nx-1).*isoold(1,2:nx-1)./cold(1,2:nx-1);		% the main event

%Northern Boundary
    c(ny,2:nx-1) = w0(ny,2:nx-1).*cold(ny,2:nx-1) +...
	wxm(ny,2:nx-1).*cold(ny,1:nx-2) +...
	wxp(ny,2:nx-1).*cold(ny,3:nx) +...
	wym(ny,2:nx-1).*cold(ny-1,2:nx-1) + consum(ny,2:nx-1);

    iso(ny,2:nx-1) = w0(ny,2:nx-1).*isoold(ny,2:nx-1) +...
	wxm(ny,2:nx-1).*isoold(ny,1:nx-2) +...
	wxp(ny,2:nx-1).*isoold(ny,3:nx) +...
	wym(ny,2:nx-1).*isoold(ny-1,2:nx-1) + Isoconsum(ny,2:nx-1).*isoold(ny,2:nx-1)./cold(ny,2:nx-1);

% The 4 corners
%southwestern boundary
    c(1,1) = w0(1,1).*cold(1,1) +...
	wxp(1,1).*cold(1,2) +...
	wyp(1,1).*cold(2,1) + consum(1,1);		% the main event

    iso(1,1) = w0(1,1).*isoold(1,1) +...
	wxp(1,1).*isoold(1,2) +...
	wyp(1,1).*isoold(2,1) + Isoconsum(1,1).*isoold(1,1)./cold(1,1);		% the main event

%southeastern boundary
    c(1,nx) = w0(1,nx).*cold(1,nx) +...
	wxm(1,nx).*cold(1,nx-1) +...
	wyp(1,nx).*cold(2,nx) + consum(1,nx);		% the main event

    iso(1,nx) = w0(1,nx).*isoold(1,nx) +...
	wxm(1,nx).*isoold(1,nx-1) +...
	wyp(1,nx).*isoold(2,nx) + Isoconsum(1,nx).*isoold(1,nx)./cold(1,nx);		% the main event

%northwestern boundary
    c(ny,1) = w0(ny,1).*cold(ny,1) +...
	wxp(ny,1).*cold(ny,2) +...
	wym(ny,1).*cold(ny-1,1) + consum(ny,1);		% the main event

    iso(ny,1) = w0(ny,1).*isoold(ny,1) +...
	wxp(ny,1).*isoold(ny,2) +...
	wym(ny,1).*isoold(ny-1,1) + Isoconsum(ny,1).*isoold(ny,1)./cold(ny,1);		% the main event

%southeastern boundary
    c(ny,nx) = w0(ny,nx).*cold(ny,nx) +...
	wxm(ny,nx).*cold(ny,nx-1) +...
	wym(ny,nx).*cold(ny-1,nx) + consum(ny,nx);		% the main event

    iso(ny,nx) = w0(ny,nx).*isoold(ny,nx) +...
	wxm(ny,nx).*isoold(ny,nx-1) +...
	wym(ny,nx).*isoold(ny-1,nx) + Isoconsum(ny,nx).*isoold(ny,nx)./cold(ny,nx);		% the main event


    Inventory(t)=sum(sum(c));
    LossTerm(t)=20;

    cold=c;					% update concentration matrix
    isoold=iso;
    
%%

while LossTerm(t)>maxError & oldyear<501
    t=t+1;
    consum=Jterm;
    Isoconsum=Jiso;
    consum(cold<2.5)=0;
    Isoconsum(cold<2.5)=0;

    cold(1:2,:)=OxyEquil*ones(2,nx);		% set oxygen boundary conditions to 250
    isoold(1:2,:)=startdel*ones(2,nx);		% set O18 boundary condition to .75

    newyear=ceil(t*dt/tyr);		% time counter
    if newyear ~= oldyear		% happy new year?
    	oldyear=newyear;		% only do once a year
        disp(oldyear); 
        figure(2)
        clabel(contour(c./250*100));
        figure(3)
        clabel(contour(((iso./c)/0.002 -1)*1000));        
        pause(.1)
    end

    c(2:ny-1,2:nx-1) = w0(2:ny-1,2:nx-1).*cold(2:ny-1,2:nx-1) +...
	wxm(2:ny-1,2:nx-1).*cold(2:ny-1,1:nx-2) +...
	wxp(2:ny-1,2:nx-1).*cold(2:ny-1,3:nx) +...
	wym(2:ny-1,2:nx-1).*cold(1:ny-2,2:nx-1) +...
	wyp(2:ny-1,2:nx-1).*cold(3:ny,2:nx-1)+ consum(2:ny-1,2:nx-1);		% the main event


    iso(2:ny-1,2:nx-1) = w0(2:ny-1,2:nx-1).*isoold(2:ny-1,2:nx-1) +...
	wxm(2:ny-1,2:nx-1).*isoold(2:ny-1,1:nx-2) +...
	wxp(2:ny-1,2:nx-1).*isoold(2:ny-1,3:nx) +...
	wym(2:ny-1,2:nx-1).*isoold(1:ny-2,2:nx-1) +...
	wyp(2:ny-1,2:nx-1).*isoold(3:ny,2:nx-1)+ Isoconsum(2:ny-1,2:nx-1).*isoold(2:ny-1,2:nx-1)./cold(2:ny-1,2:nx-1);		% the main event

%Western Boundary
    c(2:ny-1,1) = w0(2:ny-1,1).*cold(2:ny-1,1) +...
	wxp(2:ny-1,1).*cold(2:ny-1,2) +...
	wym(2:ny-1,1).*cold(1:ny-2,1) +...
	wyp(2:ny-1,1).*cold(3:ny,  1)+ consum(2:ny-1,1);		% the main event

    iso(2:ny-1,1) = w0(2:ny-1,1).*isoold(2:ny-1,1) +...
	wxp(2:ny-1,1).*isoold(2:ny-1,2) +...
	wym(2:ny-1,1).*isoold(1:ny-2,1) +...
	wyp(2:ny-1,1).*isoold(3:ny,1)+ Isoconsum(2:ny-1,1).*isoold(2:ny-1,1)./cold(2:ny-1,1);		% the main event

%Eastern Boundary
    c(2:ny-1,nx) = w0(2:ny-1,nx).*cold(2:ny-1,nx) +...
	wxm(2:ny-1,nx).*cold(2:ny-1,nx-1) +...
	wym(2:ny-1,nx).*cold(1:ny-2,nx) +...
	wyp(2:ny-1,nx).*cold(3:ny,  nx)+ consum(2:ny-1,nx);		% the main event

    iso(2:ny-1,nx) = w0(2:ny-1,nx).*isoold(2:ny-1,nx) +...
	wxm(2:ny-1,nx).*isoold(2:ny-1,nx-1) +...
	wym(2:ny-1,nx).*isoold(1:ny-2,nx) +...
	wyp(2:ny-1,nx).*isoold(3:ny,nx)+ Isoconsum(2:ny-1,nx).*isoold(2:ny-1,nx)./cold(2:ny-1,nx);		% the main event

%Southern Boundary
    c(1,2:nx-1) = w0(1,2:nx-1).*cold(1,2:nx-1) +...
	wxm(1,2:nx-1).*cold(1,1:nx-2) +...
	wxp(1,2:nx-1).*cold(1,3:nx) +...
	wyp(1,2:nx-1).*cold(2,2:nx-1)+ consum(1,2:nx-1);		% the main event

    iso(1,2:nx-1) = w0(1,2:nx-1).*isoold(1,2:nx-1) +...
	wxm(1,2:nx-1).*isoold(1,1:nx-2) +...
	wxp(1,2:nx-1).*isoold(1,3:nx) +...
	wyp(1,2:nx-1).*isoold(2,2:nx-1)+ Isoconsum(1,2:nx-1).*isoold(1,2:nx-1)./cold(1,2:nx-1);		% the main event

%Northern Boundary
    c(ny,2:nx-1) = w0(ny,2:nx-1).*cold(ny,2:nx-1) +...
	wxm(ny,2:nx-1).*cold(ny,1:nx-2) +...
	wxp(ny,2:nx-1).*cold(ny,3:nx) +...
	wym(ny,2:nx-1).*cold(ny-1,2:nx-1) + consum(ny,2:nx-1);

    iso(ny,2:nx-1) = w0(ny,2:nx-1).*isoold(ny,2:nx-1) +...
	wxm(ny,2:nx-1).*isoold(ny,1:nx-2) +...
	wxp(ny,2:nx-1).*isoold(ny,3:nx) +...
	wym(ny,2:nx-1).*isoold(ny-1,2:nx-1) + Isoconsum(ny,2:nx-1).*isoold(ny,2:nx-1)./cold(ny,2:nx-1);

% The 4 corners
%southwestern boundary
    c(1,1) = w0(1,1).*cold(1,1) +...
	wxp(1,1).*cold(1,2) +...
	wyp(1,1).*cold(2,1) + consum(1,1);		% the main event

    iso(1,1) = w0(1,1).*isoold(1,1) +...
	wxp(1,1).*isoold(1,2) +...
	wyp(1,1).*isoold(2,1) + Isoconsum(1,1).*isoold(1,1)./cold(1,1);		% the main event

%southeastern boundary
    c(1,nx) = w0(1,nx).*cold(1,nx) +...
	wxm(1,nx).*cold(1,nx-1) +...
	wyp(1,nx).*cold(2,nx) + consum(1,nx);		% the main event

    iso(1,nx) = w0(1,nx).*isoold(1,nx) +...
	wxm(1,nx).*isoold(1,nx-1) +...
	wyp(1,nx).*isoold(2,nx) + Isoconsum(1,nx).*isoold(1,nx)./cold(1,nx);		% the main event

%northwestern boundary
    c(ny,1) = w0(ny,1).*cold(ny,1) +...
	wxp(ny,1).*cold(ny,2) +...
	wym(ny,1).*cold(ny-1,1) + consum(ny,1);		% the main event

    iso(ny,1) = w0(ny,1).*isoold(ny,1) +...
	wxp(ny,1).*isoold(ny,2) +...
	wym(ny,1).*isoold(ny-1,1) + Isoconsum(ny,1).*isoold(ny,1)./cold(ny,1);		% the main event

%southeastern boundary
    c(ny,nx) = w0(ny,nx).*cold(ny,nx) +...
	wxm(ny,nx).*cold(ny,nx-1) +...
	wym(ny,nx).*cold(ny-1,nx) + consum(ny,nx);		% the main event

    iso(ny,nx) = w0(ny,nx).*isoold(ny,nx) +...
	wxm(ny,nx).*isoold(ny,nx-1) +...
	wym(ny,nx).*isoold(ny-1,nx) + Isoconsum(ny,nx).*isoold(ny,nx)./cold(ny,nx);		% the main event


    Inventory(t)=sum(sum(c));
    LossTerm(t)=Inventory(t-1)-Inventory(t);
    minOxy(newyear)=min(min(c));

    cold=c;					% update concentration matrix
    isoold=iso;
    

end
			% now must plot one final time
Oxysat=(c/250)*100;
Delfin=((iso./c)/0.002 -1)*1000;