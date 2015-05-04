function [x]=kwn_bck_sub(U,b)
n=length(U);
x=zeros(n);
x(n)=b(n)/U(n,n);
for i=n-1:-1:1
    for j=i+1:n
        b(i)=b(i)-(U(i,j)*x(j));
    end
    if U(i,i)==0
        error('Matrix U is singular')
    else
        x(i)=b(i)/U(i,i);
    end
end
x=x(:,1)';
end
