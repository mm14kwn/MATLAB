x0=0;
x1=1;
Delta=0.04;
x=x0:Delta:x1;
i=1;
y(i)=1;
while (x(i)<x1)
    y(i+1)=y(i)+(Delta)*(((y(i))^2)+1);
    i=i+1;
end
plot(x,y);