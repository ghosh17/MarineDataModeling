
function C=solarDose(t)
% Calculate photosynthesis (C) as ONLY a function of solar dose (no dependance
% on population or on nutrients
% dC/dt = q*Ed(t)
%lat = 31°40
%lon = 64°10

q=0.01;

a=0.5;
b=500;
tmp=sin(0.548)*sin(a) - cos(0.548)*cos(a).*cos((pi.*t/12)- 1.119);

tmp(tmp<0)=0;

y=b.*(0.6+0.2.*tmp).*tmp;

C=q*y;