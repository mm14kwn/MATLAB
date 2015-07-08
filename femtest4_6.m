clear all
xbpts=[0;1000;1000;0];
ybpts=[-10;-10;10;10];
poly1=[3;4;xbpts;ybpts];
d1=decsg(poly1);
[p,e,t]=initmesh(d1);
[p,e,t]=refinemesh(d1,p,e,t);
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
Nn=sum(in);
Nl=sum(left);
Nr=sum(right);
Ntop=sum(top);
Nb=sum(bottom);
NNodes=sum(in)-sum(on);
Ax=zeros(Nn);
Ay=zeros(Nn);
M=zeros(Nn);
Nt=input('Input Nt=');
Uriver=0;
U=Uriver*ones(Nn,1);
V=zeros(Nn,1);
Ulbound=Uriver*ones(Nn,1);
Urbound=zeros(Nn,1);
Vtbound=zeros(Nn,1);
Vbbound=zeros(Nn,1);
pint=NaN(2,NNodes);
pbound=NaN(2,sum(on));
pnew=NaN(2,sum(in));
tnew=NaN(size(t));
orig=NaN(NNodes,1);
orig2=NaN(sum(on),1);
l=1;
q=1;
for k=1:size(p,2)
if on(k)==0
pint(:,l)=p(:,k);
orig(l)=k;
l=l+1;
else
pbound(:,q)=p(:,k);
orig2(q)=k;
q=q+1;
end
end
pnew=[pint,pbound];
original=[orig;orig2];
for i=1:sum(in)
	originv(original(i))=i;
end
originv=transpose(originv);
for i=1:3
for j=1:size(t,2)
	tnew(i,j)=originv(t(i,j));
end
end
for i=[1,2]
	for j=1:size(e,2);
  enew(i,j)=originv(e(i,j));
end
end
tnew(4,:)=1;
enew([3,4,5,6,7],:)=e([3,4,5,6,7],:);
for i=1:max(max(t))
	h(i)=cos((pi*p(1,i))/(max(xbpts)));
	%U(i)=(2/max(xbpts))*pi*9.81*sin((2*pi*p(1,i))/(max(xbpts)));
 end
%h=zeros(max(max(t)),1);
		 tau=input('tau=');
h=transpose(h);
PSI=[U;V;h];
  psistore=zeros(length(PSI),1+Nt/tau);
psistore(:,1)=PSI;
for k=1:size(t,2)
	circ(k)=circtri(p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
for alpha=1:3
i=t(alpha,k);
for beta=1:3
j=t(beta,k);
M(i,j)=M(i,j)+mhat(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ax(i,j)=Ax(i,j)+ahatx2(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Ay(i,j)=Ay(i,j)+ahaty2(alpha,beta,p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
end

end
AAx=blkdiag(Ax,Ax,Ax);
AAy=blkdiag(Ay,Ay,Ay);
MM=blkdiag(M,M,M);
Z=zeros(length(U));
I=eye(length(U));
K1=[Z,Z,I;Z,Z,Z;I,Z,Z];
K2=[Z,Z,Z;Z,Z,I;Z,I,Z];
for time=2:1+Nt/tau
clear AAA G Gu Gv Gh
	 Gu=zeros(size(p,2),1);
Gv=zeros(size(p,2),1);
Gh=zeros(size(p,2),1);
for k=1:size(p,2)
for alpha=1:3
i=t(alpha,k);
Gu(i)=Gu(i)+Ghat(alpha,PSI(k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Gv(i)=Gv(i)+Ghat(alpha,PSI((length(PSI)/3)+k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
Gh(i)=Gh(i)+Ghat(alpha,PSI((2*length(PSI)/3)+k),p(1,t(1,k)),p(2,t(1,k)),p(1,t(2,k)),p(2,t(2,k)),p(1,t(3,k)),p(2,t(3,k)));
end
end
G=[Gu;Gv;Gh];


AAA=MM-tau*AAx*K1-tau*AAy*K2;
AAA2=AAA;
r=1;
s=1;
for i=1:length(on)
	if left(i)==1
	for j=1:length(PSI)/3
G(j)=G(j)-AAA(i,j)*Ulbound(i);
G(i)=Ulbound(i);
%G(i)=AAA(i,j)*Ulbound(i);
end
AAA(:,i)=0;
AAA(i,:)=0;
AAA(i,i)=1;
remu(r)=i;
r=r+1;
end
if right(i)==1
	for j=1:length(PSI)/3
G(j)=G(j)-AAA(i,j)*Urbound(i);
G(i)=Urbound(i);
end
AAA(i,:)=0;
AAA(:,i)=0;
AAA(i,i)=1;
remu(r)=i;
r=r+1;
end
if top(i)==1
	for j=1:length(PSI)/3
G((length(PSI)/3)+j)=G((length(PSI)/3)+j)-AAA((length(PSI)/3)+i,(length(PSI)/3)+j)*Vtbound(i);
G((length(PSI)/3)+i)=Vtbound(i);
AAA((length(PSI)/3)+i,(length(PSI)/3)+j)=0;
AAA((length(PSI)/3)+j,(length(PSI)/3)+i)=0;
end
AAA((length(PSI)/3)+i,:)=0;
AAA(:,(length(PSI)/3)+i)=0;
AAA((length(PSI)/3)+i,(length(PSI)/3)+i)=1;
remv(s)=i;
s=s+1;
end
if bottom(i)==1;
for j=1:length(PSI)/3
G((length(PSI)/3)+j)=G((length(PSI)/3)+j)-AAA((length(PSI)/3)+i,(length(PSI)/3)+j)*Vbbound(i);
G((length(PSI)/3)+i)=Vbbound(i);
AAA((length(PSI)/3)+i,(length(PSI)/3)+j)=0;
AAA((length(PSI)/3)+j,(length(PSI)/3)+i)=0;
end
AAA(:,i+length(PSI)/3)=0;
AAA(i+length(PSI)/3)=0;
AAA(i+length(PSI)/3,i+length(PSI)/3)=1;
remv(s)=i;
s=s+1;
end
end
remv=remv+(length(PSI)/3);
remove=[remu,remv];
AAA2(remove,:)=[];
AAA2(:,remove)=[];
PSIn=AAA\G;
psistore(:,time)=PSIn;
U=PSIn(1:length(PSIn)/3);
V=PSIn(1+length(PSIn)/3:2*length(PSIn)/3);
h=PSIn(1+2*length(PSIn)/3:length(PSIn));
for i=1:length(U)
	for j=1:length(U)
		E(i,j,time)=(1/2)*M(i,j)*(U(i)*U(j)+V(i)*V(j))+(1/2)*9.81*M(i,j)*(h(i)*h(j));
end
end
maxE(time)=max(max(E(:,:,time)));
minE(time)=min(min(E(:,:,time)));
PSI=PSIn;
time
end
ustore=psistore(1:length(PSIn)/3,:);
vstore=psistore(1+length(PSIn)/3:2*length(PSIn)/3,:);
hstore=psistore(1+2*length(PSIn)/3:length(PSIn),:);
%symmAAA=zeros(size(AAA2));
%for i=2:size(AAA2,1)
%	for j=1:i
%		if issymmetric(AAA2(j:i,j:i))==1
%		symmAAA(j,i)=1;
%end
%end
%end
%symmAy=zeros(size(Ay));
%for i=2:size(Ay,1)
%        for j=1:i
%                if issymmetric(Ay(j:i,j:i))==1
%                symmAy(j,i)=1;
%end
%end
%end
%symmAx=zeros(size(Ax));
%for i=2:size(Ax,1)
%        for j=1:i
%                if issymmetric(Ax(j:i,j:i))==1
%                symmAx(j,i)=1;
%end
%end
%end
