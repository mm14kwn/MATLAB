clear all
a=3600000;
b=4500000;
x=[0:1000:a];
y=[(-b/2):1000:(b/2)]; 
psi=NaN(length(y),length(x));
beta=NaN(length(y));
re=6371000;
omega=7.2921*10^(-5);
r=-0.5;
W=100;
for j=1:length(y)
    beta(j)=(2*omega/re)*cosd((10^(-5))*y(j)+37.5);
    for i=1:length(x)
        %psi(j,i)=((pi*W)/(beta(j)*b))*(cos((pi*y(j))/b))*(a*((1-exp((-beta(j)*x(i))/r))/(1-exp((-beta(j)*a)/r)))-x(i));
        psi(j,i)=((pi*W)/(beta(j)*b))*(cos((pi*y(j))/b))*(a-x(i)-a*exp(-beta(j)*x(i)/r));
    end
end
length(x)
size(psi,2) 
length(y)
size(psi,1)
contour(x,y,psi)