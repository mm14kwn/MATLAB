clear all
t0=0;
t1=40;
Delta=[0.001,0.002,0.005];
e=(10^3).*(Delta);
for j=1:3
    t=t0:Delta(j):t1;
    i=1;
    x(i)=1;
    y(i)=2;
    while (t(i)<t1)
        y(i+1)=y(i)+(Delta(j))*((y(i))*(1-2*(x(i))));
        x(i+1)=x(i)+(Delta(j))*(-(x(i))*(2-(y(i))));
        i=i+1;
    end
    plot(t,x,'b');
    hold on
    plot(t,y,'g');
    title(['{\color{blue}predator' '-' '\color{green}prey}' ' graph, \Delta=',num2str(Delta(j)),'(euler)'])
    xlabel('t')
    print('-djpeg',['(euler)_Delta',int2str(e(j))])
    hold off
    plot(x,y)
    title(['Predator vs Prey(euler), \Delta=',num2str(Delta(j))])
    xlabel('Predator')
    ylabel('Prey')
    print('-djpeg',['x_y_(euler)_Delta',int2str(e(j))])
    clear t i x y
end