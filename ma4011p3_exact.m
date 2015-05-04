function u=ma4011p3_exact(N,a,b)
h=abs(a-b)/N; %define h
x=a:h:b+h; %define x as vector of xi with spacing h and i=0:N+1
x=x'; %make vector vertical
u=-pi*x-sin(pi*x); %calculate exact solution u(xi) at all xi
end