clear all
d=[0.001 0.002 0.005];
e=d.*(10^3);
for i=1:3
    [x,y,t]=improved_euler_pp(d(i),0,40);
    plot(t,x,'b')
    hold on
    plot(t,y,'g')
    title(['{\color{blue}predator' '-' '\color{green}prey}' ' graph, \Delta=',num2str(d(i))])
    xlabel('t')
    print('-djpeg',['Delta',int2str(e(i))])
    hold off
    plot(x,y)
    title(['Predator vs Prey, \Delta=',num2str(d(i))])
    xlabel('Predator')
    ylabel('Prey')
    print('-djpeg',['x_y_Delta',int2str(e(i))])
    clear x y t
end