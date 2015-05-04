i=1;
for t=0:0.01:4
    s(i)=kwn_spline(t,[0,1,2,3,4],[0,1,0,-1,0],[0,-3,0,3,0]);
    i=i+1;
end
plot(0:0.01:4,s)
title('plot of spline values')
xlabel('t')
ylabel('s(t)')