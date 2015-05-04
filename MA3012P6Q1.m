xi=[-1,-0.5,0,0.5,1];
yi=[1.382,0.062,0.07,0.314,0.989];
c=leastsqrexp(xi,yi);
f=c(1)*exp(-1:0.01:1)+c(2)*exp(-(-1:0.01:1));
plot(xi,yi,'ob',xi,yi,'+k',-1:0.01:1,f,'r')
title(['{\color{blue}Data Points '...
'\color{red}Least Squares Approximation}'])
xlabel('x')
ylabel('y')

