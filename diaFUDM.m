function [Oxy, Iso, Oxysat, Delfin, z]= diaFUDM(depth,Oxy1, Oxy2, del1, del2, a, J, Kv, w)
%Input depth, Oxy1 -> top conc, Oxy2-> bottom O2, del1 del2 -> delO18, a->
%alpha, w-> verticle advection terms

%RECODE 1D model for dissolved oxygen consumption and isotopic
%fractionations, includes advection, diffusion and consumption
%
%writen, 12/18/2004 by Naomi Levine
%Revised 4/05 to change units to meters and seconds

%INPUT:  distance (length of model);  percentOxy (start % oxy saturation);
%del18O (start del 18 O value); a (alpha); J (consumption rate); K
%(diffusion coef.);  u (advection coeff);

%OUTPUT:  Oxy (oxygen concentration in umol/kg where satruation is 250);
%         Iso (18O concentration in umol/kg);
%         Oxysat (% oxygen saturation)
%         Delfin (del 18O values)

%UNITS:  meters and seconds.
%       J= mmol/m^3/sec
%       K= m^2/sev
%       u= m/sec
%stability:  d<=0.5, 0<=c<=1, c+2*d<=1;

yr2s=365.25*24*60*60;

dz = 50;                   %in meters... Each box
box = ceil(depth/dz); 
dz=depth/box;
z=1:dz:depth+1;%box depths 
box=length(z);

%J=1.0.*(13.1.*(z./100).^(-1.858))+0.25.*(21.3.*(z./100).^(-1.842)) +1.5.*(2.13.*(z./100).^(-1.988));
%J=-J'./yr2s;

dt= 1e6;  %same as iso           %in secs
%dt= 1e8;  %faster and still stable           %in secs

maxError= 0.00000002; %our criteria to determine steady state. 

%CONVERT Oxy and J from umol/kg to mmol/m3   7/11
J=J*1.026;

%percent saturation to concentration
startOxy= 250*Oxy1/100;
startdel= ((del1/1000+1)*0.002)*startOxy;

%percent saturation to concentration
endOxy= 250*Oxy2/100;
enddel= ((del2/1000+1)*0.002)*endOxy;

startOxy=startOxy*1.026; startdel=startdel*1.026; endOxy=endOxy*1.026; enddel=enddel*1.026;

d = Kv * dt / (dz^2);
c = w * dt / dz;
w0 = 1 - 2*d -c;
wm = d + c;
wp = d;
consum = J * dt;
Isoconsum = J*dt*a;

Oxy=ones(box, 1)*startOxy; Oxyold=Oxy;
Iso=ones(box,1)*startdel;  Isoold=Iso;
i=0;
Oxydiff=1;

while Oxydiff>maxError | i<10
%for i=1:500;
    i=i+1;
    Oxy(1) = startOxy;
    Iso(1) = enddel;
    
    Oxy(end) = endOxy;
    Iso(end) = enddel;
    
    Oxy(2:box-1) = w0 * Oxyold(2:box-1) + wm*Oxyold(1:box-2) ...
        + wp*Oxyold(3:box) + consum;
    Iso(2:box-1)= w0 * Isoold(2:box-1) + wm*Isoold(1:box-2) ...
        + wp*Isoold(3:box) + Isoconsum .* Isoold(2:box - 1)./Oxyold(2:box-1);
    
    Oxydiff=abs(sum(Oxyold-Oxy)/box);
    Oxyold=Oxy;
    Isoold=Iso;
    Oxysat=(Oxy/250)*100;
    Delfin=((Iso./Oxy)/0.002 -1)*1000;
    
    if mod(i,100)==0 
        subplot(1,2,1); plot(Oxysat, -(819.6454-z), 'k-'); subplot(1,2,2); plot(Delfin, -(819.6454-z), 'k-');
        pause(.1)
    end

end
i

%CONVERT OXY back to umol/kg
Oxy=Oxy/1.026; Iso=Iso/1.026;
Oxysat=(Oxy/250)*100;
Delfin=((Iso./Oxy)/0.002 -1)*1000;