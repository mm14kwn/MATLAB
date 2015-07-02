function out=Ghat(alpha,Unminus,x0,y0,x1,y1,x2,y2)
dc=[-1,1,0;-1,0,1];
c=[1,-1,-1;0,1,0;0,0,1];
detJ=((y2-y0)*(x1-x0))-((y0-y1)*(x0-x2));
function hout=hf(p,q)
hout=Unminus*(c(alpha,1)+c(alpha,2)*p+c(alpha,3)*q)*abs(detJ);
end
out=(1/6)*(hf(0.5,0)+hf(0,0.5)+hf(0.5,0.5));
end
