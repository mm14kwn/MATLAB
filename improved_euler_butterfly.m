function [x_out, y_out, z_out, t_out]=improved_euler_butterfly(Delta_in,T0,T)
t_max=T;
t_min=T0;
Delta=Delta_in;
t=t_min:Delta:t_max;
x=zeros(1,size(t,2));
y=zeros(1,size(t,2));
z=zeros(1,size(t,2));
w=zeros(1,3);
disp(size(t,2));
w(1)=0; 
w(2)=0;
w(3)=0;
for i=1:size(t,2)
    x(i)=w(1);
    y(i)=w(2);
    z(i)=w(3);
    w1=butterfly_rhs(w);
    w2=butterfly_rhs(w+Delta*w1);
    w=w+Delta*(w1+w2)/2;
end;
t_out=t;
x_out=x;
y_out=y;
z_out=z;

