for i=1:N
    chi(i)=(k(1)/(2*pi))*ln(y(1,i)*y(3,i)*(((x(3,i)-x(1,i))^2 + (y(1,i)+y(3,i))^2)/((x(3,i)-x(1,i))^2 + (y(3,i)-y(1,i))^2)));
end
plot(i,chi(i))
