function mat=nmat(xa,xb,ya,yb,x1,x2,x3,y1,y2,y3)
  mat=zeros(3,3);
  function ou=matint(x,y,a,b,c,d,e,f,g,h)
  ou=(1/3)*x*y*(a*e*(3*d*h-y^2)+b*f*(x^2-3*c*g)+c*g*y^2-d*h*x^2);
end
function ou2=matint2(xa2,xb2,ya2,yb2,a2,b2,c2,d2,e2,f2,g2,h2)
  ou2=matint(xa2,ya2,a2,b2,c2,d2,e2,f2,g2,h2)+matint(xb2,yb2,a2,b2,c2,d2,e2,f2,g2,h2)-matint(xa2,yb2,a2,b2,c2,d2,e2,f2,g2,h2)-matint(xb2,ya2,a2,b2,c2,d2,e2,f2,g2,h2);
end
mat(1,1)=matint2(xa,xb,ya,yb,x2,y2,x3,y3,x2,y2,x3,y3);
mat(2,2)=matint2(xa,xb,ya,yb,x1,y1,x3,y3,x1,y1,x3,y3);
mat(3,3)=matint2(xa,xb,ya,yb,x1,y1,x2,y2,x1,y1,x2,y2);
mat(1,2)=-matint2(xa,xb,ya,yb,x2,y2,x3,y3,x1,y1,x3,y3);
mat(2,1)=mat(1,2);
mat(1,3)=matint2(xa,xb,ya,yb,x2,y2,x3,y3,x1,y1,x2,y2);
mat(3,1)=mat(1,3);
mat(2,3)=-matint2(xa,xb,ya,yb,x1,y1,x3,y3,x1,y1,x2,y2);
mat(3,2)=mat(2,3);
D=det([1,x1,y1;1,x2,y2;1,x3,y3]);
mat=(1/(D^2))*mat;
end
