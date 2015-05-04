function u=ma4011p3q1(N,a,b)
h=abs(b-a)/N; %define h
e = ones(N+1,1);
A = spdiags([e -2*e e], -1:1, N+1, N+1); %create sparse matrix A
A(N+1,N+1)=-1; %set last element to -1 instead of -2
C=-(h^2)*-(pi^2)*sin(pi*(a+h:h:b)); %create vector of -h^2 * f(xi)
C=C'; %make vector vertical
c=-((h^2)/2)*-(pi^2)*sin(pi*(b)); %create final element -h^2 * f(b)/2
C=[C;c]; %add this element to the vector
u=A\C; %solve system Au=C for u
u=[0;u]; %add u(0)=0 to solution u
end