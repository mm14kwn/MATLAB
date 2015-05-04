function [z]=spline2drv(t,y)
if length(t)~=length(y)
    error('vectors are not the same length')
end
n=length(t)-1;
A=zeros(n-1,n-1);
d=NaN(1,n-1);
z=NaN(1,n-1);
for i=1:n-1
    %calculate b(i) and assign it to its place in matrix A
    A(i,i)=(t(i+2)-t(i))/3;
    %calculate d(i)
    d(i)=((y(i+2)-y(i+1))/(t(i+2)-t(i+1)))-((y(i+1)-y(i))/(t(i+1)-t(i)));
end
for i=2:n-1
    %calculate a(i) and assign it to its place in matrix A
    A(i,i-1)=(t(i+1)-t(i))/6;
end
for i=1:n-2
    %calculate c(i) and assign it to its place in matrix A
    A(i,i+1)=(t(i+2)-t(i+1))/6;
end
%Gaussian Elimination
%calculate modified values for first row
A(1,2)=A(1,2)/A(1,1);
d(1)=d(1)/A(1,1);
%apply algorithm
for i=2:n-2
    s=A(i,i)-A(i,i-1)*A(i,i+1);
    A(i,i+1)=A(i,i+1)/s;
    d(i)=(d(i)-A(i,i-1)*d(i-1))/s;
end
%calculate endpoint
d(n-1)=(d(n-1)-A(n-1,n-2)*d(n-2))/(A(n-1,n-1)-A(n-1,n-2)*A(n-2,n-1));
%back substitution
z(n-1)=d(n-1);
for i=n-2:-1:1
    z(i)=d(i)-A(i,i+1)*z(i+1);
end
z=[0,z,0];
end
    
    


