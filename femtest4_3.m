clear all
xbpts=[0;200;200;0];
ybpts=[-50;-50;50;50];
poly1=[3;4;xbpts;ybpts];
d1=decsg(poly1);
[p,e,t]=initmesh(d1);
[p,e,t]=refinemesh(d1,p,e,t);
[in,on]=inpolygon(p(1,:),p(2,:),xbpts,ybpts);
xa=min(xbpts);
xb=max(xbpts);
ya=min(ybpts);
yb=max(ybpts);
h0=10;
u0=10;
omega=4*pi;
left=zeros(size(p,2),1);
right=zeros(size(p,2),1);
bottom=zeros(size(p,2),1);
top=zeros(size(p,2),1);
for i=1:size(p,2)
	if on(i)==1
	if p(1,i)==xa
	left(i)=1;
end
if p(1,i)==xb
  right(i)=1;
end
if p(2,i)==ya
  bottom(i)=1;
end
if p(2,i)==yb
  top(i)=1;
end
end
end
Ax=zeros(size(p,2),size(p,2));
Ay=zeros(size(p,2),size(p,2));
M=zeros(size(p,2),size(p,2));
Nt=input('Input Nt=');
Uriver=0;
U=Uriver*ones(max(max(t)),1);
V=zeros(max(max(t)),1);
for i=1:max(max(t))
	h(i)=sin((pi*p(1,i))/(max(xbpts)))*cos((pi*p(2,i))/(2*max(ybpts)));
		 end
%h=zeros(max(max(t)),1);
		 tau=input('tau=');
ustore=NaN(length(U),Nt);
vstore=NaN(length(V),Nt);
hstore=NaN(length(h),Nt);
ustore(:,1)=U;
vstore(:,1)=V;
hstore(:,1)=h;
for k=1:size(t,2)
	circ(k)=circtri(p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
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
	   Gu=zeros(size(p,2),1);
Gv=zeros(size(p,2),1);
Gh=zeros(size(p,2),1);
for k=1:size(p,2)
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
%for i=1:length(h)
%	if right(i)==1
%	BBB(i)=0;
%AA(i,:)=0;
%AA(i,1)=1;
%end
%if left(i)==1
%  BBB(i)=0;
%AA(i,:)=0;
%AA(i,i)=1;
%end
%end
hn=AA\BBB;
BBu=-9.81*tau*Ax*hn+Gu;
BBv=-9.81*tau*Ay*hn+Gv;
Mu=M;
Mv=M;
%for i=1:length(U)
%	if left(i)==1
%	BBu(i)=0;
%Mu(i,:)=0;
%Mu(i,i)=1;
%end
%if right(i)==1
%  BBu(i)=0;
%Mu(i,:)=0;
%Mu(i,i)=1;
%end
%if top(i)==1
%  BBv(i)=0;
%Mv(i,:)=0;
%Mv(i,i)=1;
%end
%if bottom(i)==1
%  BBv(i)=0;
%Mv(i,:)=0;
%Mv(i,i)=1;
%end
%end
un=Mu\BBu;
vn=Mv\BBv;
ustore(:,time)=un;
vstore(:,time)=vn;
hstore(:,time)=hn;
U=un;
V=vn;
h=hn;
for i=1:length(U)
	for j=1:length(U)
		E(i,j,time)=(1/2)*M(i,j)*(U(i)*U(j)+V(i)*V(j))+(1/2)*9.81*M(i,j)*(h(i)*h(j));
end
end
maxE(time)=max(max(E(:,:,time)));
minE(time)=min(min(E(:,:,time)));
time
end
