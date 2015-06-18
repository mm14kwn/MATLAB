function [N1,N2,N3] = triN(x,y,x1,x2,x3,y1,y2,y3)
A=(1/2)*det([1,x1,y1;1,x2,y2;1,x3,y3]);
N1=(1/(2*A))*((x2*y3-x3*y2)+(y2-y3)*x+(x3-x2)*y);
N2=(1/(2*A))*((x3*y1-x1*y3)+(y3-y1)*x+(x1-x3)*y);
N3=(1/(2*A))*((x1*y2-x2*y1)+(y1-y2)*x+(x2-x1)*y);
end