clear all;
delta=0.0001;
t=0:delta:50;
x=zeros(1,length(t));
y=zeros(1,length(t));
z=zeros(1,length(t));
x(1)=0;
y(1)=1;
z(1)=0;

for i=1:(length(t)-1) 
    %I know t(i) argument is not needed in the following expressions, as
    %RKq6_dxdt and RKq6_dydt and RKq6_dzdt are not defined with t, but it means the code can
    %be used for other functions of t, x and y with minimal change.
    k_1=RKq6_dxdt(t(i),x(i),y(i),z(i));
    l_1=RKq6_dydt(t(i),x(i),y(i),z(i));
    m_1=RKq6_dzdt(t(i),x(i),y(i),z(i));
    k_2=RKq6_dxdt(t(i)+0.5*delta,x(i)+0.5*delta*k_1,y(i)+0.5*delta*l_1,z(i)+0.5*delta*m_1);
    l_2=RKq6_dydt(t(i)+0.5*delta,x(i)+0.5*delta*k_1, y(i)+0.5*delta*l_1,z(i)+0.5*delta*m_1);
    m_2=RKq6_dzdt(t(i)+0.5*delta,x(i)+0.5*delta*k_1, y(i)+0.5*delta*l_1,z(i)+0.5*delta*m_1);
    k_3=RKq6_dxdt(t(i)+0.5*delta,x(i)+0.5*delta*k_2,y(i)+0.5*delta*l_2,z(i)+0.5*delta*m_2);
    l_3=RKq6_dydt(t(i)+0.5*delta,x(i)+0.5*delta*k_2,y(i)+0.5*delta*l_2,z(i)+0.5*delta*m_2);
    m_3=RKq6_dzdt(t(i)+0.5*delta,x(i)+0.5*delta*k_2,y(i)+0.5*delta*l_2,z(i)+0.5*delta*m_2);
    k_4=RKq6_dxdt(t(i)+delta,x(i)+delta*k_3,y(i)+delta*l_3,z(i)+delta*m_3);
    l_4=RKq6_dydt(t(i)+delta,x(i)+delta*k_3,y(i)+delta*l_3,z(i)+delta*m_3);
    m_4=RKq6_dzdt(t(i)+delta,x(i)+delta*k_3,y(i)+delta*l_3,z(i)+delta*m_3);
    x(i+1)=x(i)+((delta*(k_1+2*k_2+2*k_3+k_4))/6);
    y(i+1)=y(i)+((delta*(l_1+2*l_2+2*l_3+l_4))/6);
    z(i+1)=z(i)+((delta*(m_1+2*m_2+2*m_3+m_4))/6);
end

u=zeros(1,length(t));
v=zeros(1,length(t));
w=zeros(1,length(t));
u(1)=0.001;
v(1)=1;
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

disp(['Numerical solution: x(0)=0 (x,y,z) = ',mat2str([x(length(t)-1),y(length(t)-1),z(length(t)-1)])])
disp(['Numerical solution: x(0)=0.001 (x,y,z) = ',mat2str([u(length(t)-1),v(length(t)-1),w(length(t)-1)])])
plot3(x,y,z,'b')
grid on
hold on
plot3(u,v,w,'r')
xlabel('x')
ylabel('y')
zlabel('z')
title({'Runge Kutta approximation to butterfly effect, \Delta=0.001';'\color{red}x(0)=0.001, \color{blue}x(0)=0'})
hold off
print('-djpeg','DEAD_P6_Q4')
