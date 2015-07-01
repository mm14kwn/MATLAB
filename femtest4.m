clear all
poly1=[2;10;0;0;5;15;20;25;25;20;15;5;10;-10;-5;-5;-1;-1;1;1;5;5];
d1=decsg(poly1);
[p,e,t]=initmesh(d1);
Ax=zeros(size(t,2),size(t,2));
Ay=zeros(size(t,2),size(t,2));
for k=1:size(t,2)
for alpha=1:3
i=t(alpha,k);
for beta=1:3
j=t(beta,k);
Ax(i,j)=Ax(i,j)+Ahatx(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ay(i,j)=Ay(i,j)+Ahaty(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
