function z=Stommel2(x,y)
psi0=10^9; % setting the stream value scale
lambda=10^9;
b=2*pi*10^8;
d=20000;
gamma=1;
fprime=2*10^(-13);
R=0.02;
A=(gamma*pi*lambda^2)/(R*b*psi0);
alpha=lambda*fprime*d/R;
mu=pi*lambda/b;
m1=(-alpha-sqrt(alpha^2+4*mu^2))/2;
m2=(-alpha+sqrt(alpha^2+4*mu^2))/2;
c1=(exp(m2)-1)/(exp(m2)-exp(m1));
c2=(1-exp(m1))/(exp(m2)-exp(m1));
[x,y]=meshgrid(0:0.01:1,0:b/lambda/100:b/lambda);
z=A*b^2/(lambda^2*pi^2)*(-1+c1*exp(m1*x)+c2*exp(m2*x)).*sin(pi*lambda*y/b);
contour(x,y,z)