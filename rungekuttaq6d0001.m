clear all;
delta=0.001;
t=0:delta:400;
u=zeros(1,length(t));
v=zeros(1,length(t));
w=zeros(1,length(t));
u(1)=0.001;
v(1)=0;
w(1)=0;

for i=1:(length(t)-1) 
    %I know t(i) argument is not needed in the following expressions, as
    %RKq6_dxdt and RKq6_dydt and RKq6_dzdt are not defined with t, but it means the code can
    %be used for other functions of t, x and y with minimal change.
    p_1=RKq6_dxdt(t(i),u(i),v(i),w(i));
    q_1=RKq6_dydt(t(i),u(i),v(i),w(i));
    r_1=RKq6_dzdt(t(i),u(i),v(i),w(i));
    p_2=RKq6_dxdt(t(i)+0.5*delta,u(i)+0.5*delta*p_1,v(i)+0.5*delta*q_1,w(i)+0.5*delta*r_1);
    q_2=RKq6_dydt(t(i)+0.5*delta,u(i)+0.5*delta*p_1, v(i)+0.5*delta*q_1,w(i)+0.5*delta*r_1);
    r_2=RKq6_dzdt(t(i)+0.5*delta,u(i)+0.5*delta*p_1, v(i)+0.5*delta*q_1,w(i)+0.5*delta*r_1);
    p_3=RKq6_dxdt(t(i)+0.5*delta,u(i)+0.5*delta*p_2,v(i)+0.5*delta*q_2,w(i)+0.5*delta*r_2);
    q_3=RKq6_dydt(t(i)+0.5*delta,u(i)+0.5*delta*p_2,v(i)+0.5*delta*q_2,w(i)+0.5*delta*r_2);
    r_3=RKq6_dzdt(t(i)+0.5*delta,u(i)+0.5*delta*p_2,v(i)+0.5*delta*q_2,w(i)+0.5*delta*r_2);
    p_4=RKq6_dxdt(t(i)+delta,u(i)+delta*p_3,v(i)+delta*q_3,w(i)+delta*r_3);
    q_4=RKq6_dydt(t(i)+delta,u(i)+delta*p_3,v(i)+delta*q_3,w(i)+delta*r_3);
    r_4=RKq6_dzdt(t(i)+delta,u(i)+delta*p_3,v(i)+delta*q_3,w(i)+delta*r_3);
    u(i+1)=u(i)+((delta*(p_1+2*p_2+2*p_3+p_4))/6);
    v(i+1)=v(i)+((delta*(q_1+2*q_2+2*q_3+q_4))/6);
    w(i+1)=w(i)+((delta*(r_1+2*r_2+2*r_3+r_4))/6);
end
plot3(u,v,w)
grid on
disp(['Numerical solution: (u,v,w) = ',mat2str([u(length(t)-1),v(length(t)-1),w(length(t)-1)])])