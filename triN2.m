function [N1,N2,N3] = triN2(x,y,x1,x2,x3,y1,y2,y3)
Nh=zeros(3);
xm=[x1,x2,x3];
ym=[y1,y2,y3];
for j=1:3
    kmat=1:3;
    kmat(j)=[];
    for k=kmat
        lmat=1:3;
        ind=[j,k];
        lmat(ind)=[];
        for l=lmat
            D=det([1,x,y;1,xm(k),ym(k);1,xm(l),ym(l)]);
            C=det([1,xm(j),ym(j);1,xm(k),ym(k);1,xm(l),ym(l)]);
            Nh(j)=Nh(j)+D/C;
        end
    end
end
N1=Nh(1)/2;
N2=Nh(2)/2;
N3=Nh(3)/2;
end