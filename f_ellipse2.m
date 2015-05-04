function z=f_ellipse2(x,y)
a=1;
b=2;
c=3;
z=a.*b.*c.*cos(x).*sin(y).*sin(x).*sin(y).*cos(y).*sin(y).*(sqrt(((((b.^2).*(c.^2).*((cos(x)).^2))+((a.^2).*(c.^2).*((sin(x)).^2))).*((sin(y)).^2))+((a.^2).*(b.^2).*((cos(y)).^2))));
end