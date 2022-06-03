function [a,sa,chisqr,covmat] = surfit(x,y,z,n)
%
%		[a, sa, chisqr, covmat] = surfit(x,y,z,n)
%
%	returns the parameters of a nth order fit (a in order
%	of const, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3 etc)
%	for values of n=0 on up, along with  
%	errors (sa), chi squared and the covariance matrix for
%	the data stored in x, y,  and z (N-length 
%	column vectors).  Routine uses SVD to solve normal equations.
%
M=.5*n*(n+3)+1;				% this is number of parameters
%		first check out inputs
%			x.....

[N L]=size(x);
if L ~= 1, x=x'; [N L]=size(x); end
if N < M, error('Insufficient number of data for fit'); end

%			and then y.....

[Ny L]=size(y);
if L ~= 1, y=y'; [Ny L]=size(y); end
if N ~= Ny, error('x and y must be same length'); end

%			and finally z...

[Nz L]=size(z);
if L ~= 1, z=z'; [Nz L]=size(z); end
if N ~= Nz, error('x and z must be same length'); end
%
%	 Now construct design matrix
%
A=[ones(N,1)];
if n>0
	for i=1:n
	for j=0:i
		A = [A (x.^(i-j)).*(y.^j) ];
	end
	end
end
%
%		then do fit
%
[U S V]=svd(A,0);		% Then solve by SVD
w=diag(S);			% Get singular values
wmin=max(w)*1.e-5;		% find minimum s.v. permissable
for i=1:M			% and invert diagonal
	if w(i) > wmin
		w(i)=1/w(i);
	    else
		w(i)=0;
	end
end
Si=diag(w);			% and remake matrix
a=V*Si*(U'*z);			% compute coefficient vector
covmat=V*Si^2*V';		% compute covariance matrix
chisqr=sum((A*a-z).^2)/(N-M);	% calculate chisquare
sa=sqrt(chisqr*diag(covmat));	% errors from chisqr * diag of cov matrix
