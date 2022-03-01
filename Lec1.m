A= [3 5 5; 2 4 2; 3 10 10];

B = [10 3 3 ; 5 5 7; 2 6 8];

A*B;

B*A;

A.*B;

A = [1 2 1; 2 3 -2; 1 -2 3];

b = [8; 2; 6];

x1 = inv(A)*b;

x2 = A\b;

eye(3); %identity matrix of 3X3

v1= [4 3];
v2 = [5 7];

%%%%%%%%%%%%%%
%Rank Deficient 

A = [1 2 1; 2 3 -2; 3 5 -1];

b = [8; 2; 10];

x = A\b;

%rank(A); %This rank is 2 becasue it is rank deficient coz third equation is a sum of first two and is therefore not linearly independent

%det = det(A);

[U,S,V] = svd(A,0)

%U*S*V'

%Now find out inv(A) or BGI

s = diag(S)

%size of vector s

n = size(s);

num = n(:,1);

%find 0 element in s vector
itr = 1;
position = -1;
while (itr<=num)
    if s(itr) <= 1e-10
        s(itr) = 1;
        position = itr;
    end
    itr = itr + 1;
end

w = 1./s;

if(position > -1)
    w(position) = 0; 
end

W = diag(w);

BGI = V*W*U'

x = BGI*b

A*BGI %notice how it is not I

cond(A) %if cond is very large then matrix is more singular. FUction RCOND is the reciprocal of this fucntion

