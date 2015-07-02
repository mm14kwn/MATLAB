function out=Ahatx(alpha,beta,x0,y0,x1,y1,x2,y2)
d=[-1,1,0;-1,0,1];
c=[1,-1,-1;0,1,0;0,0,1];
detJ=((y2-y0)*(x1-x0))-((y0-y1)*(x0-x2));
function hout=hf(p,q)
hout=(c(alpha,1)+c(alpha,2)*p+c(alpha,3)*q)*((y2-y0)*d(1,beta)+(y0-y2)*d(2,beta))*(1/detJ)*abs(detJ);
end
out=(1/6)*(hf(0.5,0)+hf(0,0.5)+hf(0.5,0.5));
end
