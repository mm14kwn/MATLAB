function[gi]=vortexg(i,k,x,y,N)
gi=0;
jmat=[1:1:N];
jmat(i)=[];
for j=jmat
    rr=((x(i)-x(j))^2)+((y(i)-y(j))^2);
    gi=gi+(k(j)*(x(i)-x(j))/rr);
end
end