function [Arot] = varimax(Ar)
%   Arot = varimax(Ar)
%       performs a varimax rotation on the factor loading
%       matrix "Ar", returning the rotated loading matrix "Arot"
% This procedure follows algorithm as spelled out in
% Harman (1960) in Chapter 14, section 4. 
MaxIts = 25;        % maximum number of iterations
Tolerance = 0.0001;     % convergence criterion
b=Ar;
[n,nf]=size(Ar);
hjsq=diag(Ar*Ar');   % communalities
hj=sqrt(hjsq);
bh=Ar./(hj*ones(1,nf));     % compute variances of loadings
V0=n*sum(sum(bh.^4))-sum(sum(bh.^2).^2);
for it=1:MaxIts; % Maximum number of iterations
for i=1:nf-1 % Program cycles through 2 factors
  jl=i+1;    % at a time.
  for j=jl:nf
      xj=Ar(:,i)./hj;   % notation here closely
      yj=Ar(:,j)./hj;   % follows harman
      uj=xj.*xj-yj.*yj;
      vj=2*xj.*yj;
      A=sum(uj);
      B=sum(vj);
      C=uj'*uj-vj'*vj;
      D=2*uj'*vj;
      num=D-2*A*B/n;
      den=C-(A^2-B^2)/n;
      tan4p=num/den;
      phi=atan2(num,den)/4;
      angle=phi*180/pi;
      [i j it angle];
      if abs(phi)>.00001;
          Xj=cos(phi)*xj+sin(phi)*yj;
          Yj=-sin(phi)*xj+cos(phi)*yj;
          bj1=Xj.*hj;
          bj2=Yj.*hj;
          b(:,i)=bj1;
          b(:,j)=bj2;
          Ar(:,i)=b(:,i);
          Ar(:,j)=b(:,j);
      end
  end
end;
Arot=b;
bh=Ar./(hj*ones(1,nf));     % compute variances of loadings
V=n*sum(sum(bh.^4))-sum(sum(bh.^2).^2);
if abs(V-V0)<Tolerance;break;else V0=V;end; % converged?
end;

