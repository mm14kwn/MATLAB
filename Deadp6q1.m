clear all
r0=0;
r1=400;
Delta=0.001;
a=0.2;
b=0.2;
c=5.7;
r=r0:Delta:r1;
i=1;
x(i)=0;
y(i)=0;
z(i)=0;
while (r(i)<r1)
    x(i+1)=x(i)+(Delta)*(-y(i)-z(i));
    y(i+1)=y(i)+(Delta)*(x(i)+a*y(i));
    z(i+1)=z(i)+(Delta)*(b+(z(i))*(x(i)-c));
    i=i+1;
end
plot3(x,y,z)
grid on