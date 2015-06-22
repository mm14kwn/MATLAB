clear all
poly1=[2;10;0;0;5;15;20;25;25;20;15;5;10;-10;-5;-5;-1;-1;1;1;5;5];
d1=decsg(poly1);
[p,e,t]=initmesh(d1);
hx=0.1;
hy=0.2;
x=0:hx:25;
y=-10:hy:10;
xy=[x',y']
Nx=length(x);
Ny=length(y);
for i=1:size(p,2)
    f(i)=func(p(1,i),p(2,i));
end
for k=1:size(t,2)
for l=1:3)
trix(l)=p(1,t(l,k));
triy(l)=p(2,t(l,k));
end
in(k)=inpolygon(x,y,trix,triy);
end
%for i=1:Nx
    %for j=1:Ny
