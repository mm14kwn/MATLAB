function out=Ahatx2(alpha,beta,x0,y0,x1,y1,x2,y2)
d=[-1,1,0;-1,0,1];
c=[1,-1,-1;0,1,0;0,0,1];
detJ=((y2-y0)*(x1-x0))-((y1-y0)*(x2-x0));
function hout=hf(p,q)
lam=c(alpha,1)+c(alpha,2)*p+c(alpha,3)*q;
lambda=[lam;0];
JT=[x1-x0,y1-y0;x2-x0,y2-y0];
D=[d(1,beta);d(2,beta)];
hout=lambda.*(inv(JT)*D)*abs(detJ);
end
out=(1/6)*(hf(0.5,0)+hf(0,0.5)+hf(0.5,0.5));
end
