clear all
pnum=input('plot number = ');
N=input('N=x range = ');
M=input('M=y range = ');
dx=1/N;
dy=dx;
dt=0.001;
tau=input('tau = ');
fprime=input('fprime = ');
tidal_period=input('tidal period = ');
a=zeros(1,M);
u=zeros(M+1,N+1);
v=zeros(M+1,N+1);
eta=zeros(M,N);
a1=input('lower edge of inlet 1 (proportion of M) = ');
a2=input('upper edge of inlet 1 (proportion of M) = ');
a3=input('lower edge of inlet 2 (proportion of M) = ');
a4=input('upper edge of inlet 2 (proportion of M) = ');
phase1=input('phase shift of inlet 1 = ');
phase2=input('phase shift of inlet 2 = ');
amp=input('amplitude of tide = ');
a(a1*M:a2*M)=amp;
a(a3*M:a4*M)=amp;
t=0;
T=input('max time? ');
while t<T
  for j=1:M
  for i=2:N
  px(j,i)=(eta(j,i)-eta(j,i-1))/dx;
  end
  for i=2:N
  vf(j,i)=(v(j,i)+v(j+1,i)+v(j,i-1)+v(j+1,i-1))/4;
  u(j,i)=u(j,i)+dt*(-px(j,i)+fprime*vf(j,i)-u(j,i)/tau);
  end
  end
  for j=a1*M:a2*M
  px(j,1)=(eta(j,1)-a(j)*sin((2*pi*t/tidal_period)+phase1))/dx;
  end
  for j=a3*M:a4*M
  px(j,1)=(eta(j,1)-a(j)*sin((2*pi*t/tidal_period)+phase2))/dx;
  end
  for j=1:M
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
t=t+dt;
end
fileID=fopen(['E:\Google Drive\Lecture Notes\Leeds Semester 2\MATH5458 Geophysical Fluids\Project\parameters_',int2str(pnum),'.txt'],'w');
fprintf(fileID,'N=%4g \n',N);
fprintf(fileID,'M=%4g \n',M);
fprintf(fileID,'tau=%4g \n',tau);
fprintf(fileID,'fprime=%4g \n',fprime);
fprintf(fileID,'tidal period=%4g\n',tidal_period);
fprintf(fileID,'1st inlet=[%4g,%4g]*M \n',a1,a2);
fprintf(fileID,'2nd inlet=[%4g,%4g]*M \n',a3,a4);
fprintf(fileID,'phase shift of 1st inlet=%4g \n',phase1);
fprintf(fileID,'phase shift of 2nd inlet=%4g \n',phase2);
fprintf(fileID,'amplitude=%4g \n',amp);
fprintf(fileID,'final time=%4g \n',T);
fclose(fileID);
[C,h]=contour(eta,'k');
hold on
clabel(C,h);
quiver(u,v,'k');
hold off
print(['E:\Google Drive\Lecture Notes\Leeds Semester 2\MATH5458 Geophysical Fluids\Project\contour_',int2str(pnum),'.eps'],'-depsc')
clf
surf(eta);
print(['E:\Google Drive\Lecture Notes\Leeds Semester 2\MATH5458 Geophysical Fluids\Project\surf_',int2str(pnum),'.eps'],'-depsc')