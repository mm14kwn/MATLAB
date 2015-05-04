h=input('Input size increment: ');
L=input('Input x/y length: ');
D=input('Input depth: ');
U=input('Input U: ');
Y=[0:h:L];
Z=[0:-h:-D];
size(Y)
size(Z)
u=zeros(L/h,D/h);
for i=1:(L/h)+1
    for j=1:(D/h)+1
        u(i,j)=U*sin(pi*Y(i)/L)*cos((pi*Z(j))/(2*D));
    end
end
size(u)
contour(Y,Z,u')