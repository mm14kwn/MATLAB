clear all
clf
N=input('N=');
dx=1/N;
dt=0.1*dx;
amp=input('amplitude=');
eta=zeros(1,N);
u=zeros(1,N+1);
u(N+1)=-1;
tau=input('linear damping time=');
tidal_period=input('tidal period =');
t=0;
T=input('max time?');
etat=NaN(T/dt,N);
k=1;
while t<T
for i=2:N-1;
    px(i)=(eta(i)-eta(i-1))/dx;
end
px(1)=(eta(1)-amp*sin(2*pi*t/tidal_period))/dx;
px(N)=(1-eta(N))/dx;
u(1:N)=u(1:N)+dt*(-px-(u(1:N)/tau));
for i=1:N
    ux(i)=(u(i+1)-u(i))/dx;
end
eta=eta-(dt*ux);
t=t+dt
k=k+1;
etat(k,:)=eta;
end
for i=1:k
    plot(1:N,etat(i,:));
    axis([1 N -1.5 1.5])
    i*dt
    pause(0.001)
end
    