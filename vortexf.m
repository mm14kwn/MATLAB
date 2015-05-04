function[fi]=vortexf(i,k,x,y,N)
fi=0;
jmat=[1:1:N];
jmat(i)=[];
for j=jmat
    rr=((x(i)-x(j))^2)+((y(i)-y(j))^2);
    fi=fi-(k(j)*(y(i)-y(j))/rr);
end
end