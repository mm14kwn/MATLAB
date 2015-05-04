function c = leastsqrexp(xi, yi) 
if length(xi)~=length(yi)
    error('please input xi,yi the same length')
end
n=length(xi);
D=[0,n;n,0];
e(1,1)=dot(yi,exp(xi));
e(2,1)=dot(yi,exp(-xi));
for k=1:n
    D(1,1)=D(1,1)+exp(2*xi(k));
    D(2,2)=D(2,2)+exp(-2*xi(k));
end
c=D\e;
end
