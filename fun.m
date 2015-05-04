function[f]=fun(t,k,add,pos)
fx=zeros(4,1);
fy=zeros(4,1);
for i=1:4
    for j=1:4
        if i~=j
            fx(i)=fx(i)-Funx(t,k,i,j,add,pos);
            fy(i)=fy(i)+Funy(t,k,i,j,add,pos);
        end
    end
end
f=[fx,fy];