function w_out=butterfly_rhs(w)
p=zeros(1,3);
a=0.2;
b=0.2;
c=5.7;
p(1)=-w(2)-w(3);
p(2)=(w(1))+a*(w(2));
p(3)=b+(w(3))*(w(1)-c);
w_out=p;
