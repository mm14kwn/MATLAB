clear all
gridspacex=10;
gridspacey=10;
hx=1;
hy=1;
Nx=400;
Ny=100;
xi=-Nx/2:gridspacex:Nx/2;
yi=-Ny/2:gridspacey:Ny/2;
x=-Nx/2:hx:Nx/2;
y=-Ny/2:hy:Ny/2;
u=NaN(length(x),length(y));
err=NaN(length(x),length(y));
f=NaN(length(x),length(y));
perr=NaN(length(x),length(y));
nodes=cell(length(x),length(y));
for i=1:1+Nx/hx
    for j=1:1+Ny/hy
        for k=1:Nx/gridspacex
            if (x(i)<xi(k+1))&&(x(i)>=xi(k))
                xnode=k;
                break
            end
        end
        for l=1:Ny/gridspacey
            if y(j)<(yi(l+1))&&(y(j)>=yi(l))
                ynode=l;
                break
            end
        end
        if (x(i)/gridspacex)>=(y(j)/gridspacey)
            uplow=1;
        else
            uplow=0;
        end
        x1=xi(xnode);
        x2=xi(xnode+1);
        x3=xi(xnode+uplow);
        y1=yi(ynode);
        y2=yi(ynode+1-uplow);
        y3=yi(ynode+1);
        [Na,Nb,Nc]=triN2(x(i),y(j),x1,x2,x3,y1,y2,y3);
        u(i,j)=Na*func(x1,y1)+Nb*func(x2,y2)+Nc*func(x3,y3);
        f(i,j)=func(x(i),y(j));
        err(i,j)=u(i,j)-func(x(i),y(j));
        perr(i,j)=err(i,j)/func(x(i),y(j));
        nodes{i,j}=[x(i),y(j);x1,y1;x2,y2;x3,y3];
    end
end
contour(u)
print('-dpdf','femtestu.pdf')
contour(f)
print('-dpdf','femtestf.pdf')
exit