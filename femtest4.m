clear all
poly1=[3;4;0;100;100;0;0;0;10;10];
d1=decsg(poly1);
[p,e,t]=initmesh(d1,'hmax',50);
Ax=zeros(size(p,2),size(p,2));
Ay=zeros(size(p,2),size(p,2));
M=zeros(size(p,2),size(p,2));
Nt=10;
Uriver=1;
U=Uriver*ones(size(p,2),1);
V=zeros(size(p,2),1);
h=zeros(size(p,2),1);
tau=1;
ustore=NaN(size(p,2),Nt);
vstore=NaN(size(p,2),Nt);
hstore=NaN(size(p,2),Nt);
ustore(:,1)=U;
vstore(:,1)=V;
hstore(:,1)=h;
for k=1:size(t,2)
for alpha=1:3
i=t(alpha,k);
for beta=1:3
j=t(beta,k);
i
j
M(i,j)=M(i,j)+mhat(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ax(i,j)=Ax(i,j)+ahatx(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ay(i,j)=Ay(i,j)+ahaty(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
end
end
Minv=inv(M);
for time=2:Nt
	   Gu=zeros(size(p,2),1);
Gv=zeros(size(p,2),1);
Gh=zeros(size(p,2),1);
for k=1:size(t,2)
for alpha=1:3
i=t(alpha,k);
Gu(i)=Gu(i)+Ghat(alpha,U(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Gv(i)=Gv(i)+Ghat(alpha,V(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Gh(i)=Gh(i)+Ghat(alpha,h(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
end
AA=M-(tau^2)*(Ax*Minv*Ax+Ay*Minv*Ay);
BB=-tau*(Ax*Minv*Gu-Ay*Minv*Gv);
BBB=Gh-BB;
hn=AA\BBB
BBu=-tau*Ax*Minv*Ax*hn+Ax*Minv*Gu;
BBv=-tau*Ay*Minv*Ay*hn+Ay*Minv*Gv;
un=Ax\BBu;
vn=Ay\BBv;
ustore(:,time)=un;
vstore(:,time)=vn;
hstore(:,time)=hn;
U=un;
V=vn;
h=hn;
end
