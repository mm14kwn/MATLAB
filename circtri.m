function r=circtri(x1,y1,x2,y2,x3,y3)
xa = x2-x1;
ya = y2-y1;
xb = x3-x1;
yb = y3-y1;
da = xa^2+ya^2;
db = xb^2+yb^2;
A = 2*(xa*yb-ya*xb);
r = sqrt(da*db*((x3-x2)^2+(y3-y2)^2))/abs(A);
end
