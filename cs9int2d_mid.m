clear all;                                              %    clears all variables 
format long;                                         %    outputs numbers in 15 decimal digits 
m = 100;                                               %     number of subintervals in the x-direction 
n = 50;                                                %     number of subintervals in the y-direction 
a = 0; b = 2*pi; c = 0; d = pi;                    %     define the endpoints of integration 
hx = (d-c)/m; hy = (b-a)/n;                   %     set up the lengths hx and hy 
x = linspace(c,d,m+1);                        %     create the vector (x1; x2; : : : ; xm+1) 
y = linspace(a,b,n+1);                          %    create the vector (y1; y2; : : : ; yn+1) 
a=0;
for i=1:m
    for j=1:n                                     %calculate the approximation to the integral, 
      a=a+f_ellipse((x(i)+x(i+1))/2,(y(j)+y(j+1))/2,1,2,3);
    end
end
integral=hx*hy*a
