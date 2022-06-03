%DiaDriver:  to run diaFUDM_fn for a range of alpha and J starting for 3
%longitudinal groups
%Written 4/06 
%innitially use pre-definied w (vertical advection) and Kv (vertical
%diffusion)
%UNITS:  meters and seconds.
%       J= umol/kg/sec
%       Kv= m^2/sev
%       uv= m/sec
%OXY AND J ARE BOTH IN umol/kg CONVERTED TO mmol/m3 IN THE MODEL FN.

% addpath /Users/naomimarcil/Mfiles/seawater
load Data_sta

yr2s=365.25*24*60*60;
a=[0.981]; %0.990
J=[-.4]./yr2s;          %for u=0  K=e-5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eastern group 1.9W to 10.4E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Group3a=Group3(Group3(:,w_press)>250&Group3(:,w_press)<900, :);
Group3b=Group3(Group3(:,w_press)>=940, :);

DATAOxy=Group3a(:,w_winkler)./Group3a(:,w_oxysatKG)*100;
DATADel=Group3a(:,w_del);
DATAdepth=sw_dpth(Group3a(:,w_press),nanmean(Group3a(:,w_lat)));

[Depthsort, I]=sort(DATAdepth, 1);
Delsort=DATADel(I);
Oxysort=DATAOxy(I);

Oxy2=Oxysort(1);
del2=Delsort(1);
Oxy1=Oxysort(end);
del1=Delsort(end);

lat2=nanmean(Group3a(:,w_lat));
Dep2=Depthsort(1);
Dep1=Depthsort(end);
depth=Dep1-Dep2;

%Set w and Kv values 
w=1e-8;  %1e-8;
Kv=1E-5;

figure(1)
subplot(1,2,1); plot(Oxysort, -Depthsort, 'b*'); xlabel('%Oxy sat'); ylabel('Depth (m)'); hold on
subplot(1,2,2); plot(Delsort, -Depthsort, 'b*'); xlabel('\delta18O'); ylabel('Depth (m)'); hold on
[Oxy, Iso, Oxysat, Delfin, z]= diaFUDM(depth,Oxy1, Oxy2, del1, del2, a, J, Kv, w); 
ModDep=(Dep1-z);

figure(1)
subplot(1,2,1); plot(Oxysat, -ModDep, 'r-','linewidth', 2); xlabel('%Oxy sat'); ylabel('Depth (m)')
subplot(1,2,2); plot(Delfin, -ModDep, 'r-','linewidth', 2); xlabel('\delta18O'); ylabel('Depth (m)')

figure(1)
subplot(1,2,1); plot(Oxysat, -ModDep, 'r-', Oxysort, -Depthsort, 'b*'); xlabel('%Oxy sat'); ylabel('Depth (m)'); hold on
subplot(1,2,2); plot(Delfin, -ModDep, 'r-', Delsort, -Depthsort, 'b*'); xlabel('\delta18O'); ylabel('Depth (m)');hold on