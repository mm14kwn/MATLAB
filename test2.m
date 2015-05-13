N=100;
M=50;
dx=1/N;
dy=dx;
dt=0.001;
tau=20;
fprime=pi;
tidal_period=5;
a=zeros(1,M);
u=zeros(M+1,N+1);
v=zeros(M+1,N+1);
eta=zeros(M+1,N+1);
a(1:M/5)=0.05;
a(4*M/5:M)=-0.05;
t=0;
while t<2.5
  for j=1:M
  for i=2:N
  px(j,i)=(eta(j,i)-eta(j,i-1))/dx;
  end
  for i=2:N
  vf(j,i)=(v(j,i)+v(j+1,i)+v(j,i-1)+v(j+1,i-1))/4;
  u(j,i)=u(j,i)+dt*(-px(j,i)+fprime*vf(j,i)-u(j,i)/tau);
  end
  end
  for j=1:M
  px(j,1)=(eta(j,1)-a(j)*sin(2*pi*t/tidal_period))/dx;
  vf(j,1)=(v(j,1)+v(j+1,1))/2;
  u(j,1)=u(j,1)+dt*(-px(j,1)-u(j,1)/tau);
  end
  for i=2:N
  for j=2:M
  py(j,i)=(eta(j,i)-eta(j-1,i))/dy;
  uf(j,i)=(u(j,i)+u(j,i+1)+u(j-1,i)+u(j-1,i+1))/4;
v(j,i)=v(j,i)+dt*(-py(j,i)-fprime*uf(j,i)-v(j,i)/tau);
  end
  end
  for j=1:M
  for i=1:N
	  eta(j,i)=eta(j,i)-dt*((u(j,i+1)-u(j,i))/dx+(v(j+1,i)-v(j,i))/dy);
end
end
u(M+1,:)=u(M,:);
u(M+1,:)=0;
v(:,N+1)=v(:,N);
v(:,N+1)=0;
t=t+dt
end
