function [x_out, y_out, t_out]=improved_euler_pp(Delta_in,T0,T)
t_max=T;
t_min=T0;
Delta=Delta_in;
t=t_min:Delta:t_max;
x=zeros(1,size(t,2));
y=zeros(1,size(t,2));
z=zeros(1,2);
z(1)=1; 
z(2)=2;
for i=1:size(t,2)
    x(i)=z(1);
    y(i)=z(2);
    z1=pred_prey_rhs(z);
    z2=pred_prey_rhs(z+Delta*z1);
    z=z+Delta*(z1+z2)/2;
end;
t_out=t;
x_out=x;
y_out=y;