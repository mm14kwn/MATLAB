clear all
poly1=[2;10;0;0;5;15;20;25;25;20;15;5;10;-10;-5;-5;-1;-1;1;1;5;5];
d1=decsg(poly1);
[p,e,t]=initmesh(d1);
Ax=zeros(size(t,2),size(t,2));
Ay=zeros(size(t,2),size(t,2));
M=zeros(size(t,2),size(t,2));
Nt=100;
Uriver=1;
U=Uriver*ones(1,size(t,2));
V=zeros(1,size(t,2));
h=zeros(1,size(t,2));
tau=1;
ustore=NaN(Nt,size(t,2));
vstore=NaN(Nt,size(t,2));
hstore=NaN(Nt,size(t,2));
ustore(1,:)=U;
vstore(1,:)=V;
hstore(1,:)=h;
for k=1:size(t,2)
for alpha=1:3
i=t(alpha,k);
for beta=1:3
j=t(beta,k);
M(i,j)=M(i,j)+mhat(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ax(i,j)=Ax(i,j)+ahatx(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ay(i,j)=Ay(i,j)+ahaty(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
end
end
Minv=inv(M);
for time=2:Nt
Gu=zeros(1,size(t,2));
Gv=zeros(1,size(t,2));
Gh=zeros(1,size(t,2));
for k=1:size(t,2)
for alpha=1:3
i=t(alpha,k);
Gu(i)=Gu(i)+Ghat(alpha,U(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Gv(i)=Gv(i)+Ghat(alpha,V(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Gh(i)=Gh(i)+Ghat(alpha,H(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
end
AA=M-(tau^2)*(Ax*Minv*Ax+Ay*Minv*Ay);
BB=-tau*(Ax*Minv*Gu-Ay*Minv*Gv);
BBB=Gh-BB;
hn=AA\BBB;
BBu=-tau*Ax*Minv*Ax*hn+Ax*Minv*Gu;
BBv=-tau*Ay*Minv*Ay*hn+Ay*Minv*Gv;
un=Ax\BBu;
vn=Ay\BBv;
ustore(time,:)=un;
vstore(time,:)=vn;
hstore(time,:)=hn;
U=un;
V=vn;
h=hn;
end
