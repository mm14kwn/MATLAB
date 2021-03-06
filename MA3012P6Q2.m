clear all
k=1;
for h=[0.0025,0.005,0.01,0.02,0.04,0.08]
    %set initial conditions
    x(1)=0.1;
    ge(1)=0;
    fu(1)=10/9;
    dfu(1)=100/81;
    i=2;
    %calculate values at each point
    for t=h:h:0.4
        x(i)=x(i-1)+h*(1/(1-x(i-1)));
        fu(i)=1/(1-x(i));
        dfu(i)=1/((x(i)-1)^2);
        ge(i)=abs((1-sqrt(0.81-2*t))-x(i));
        i=i+1;
    end
    %calculate local error
    le=abs((h^2)*(1/2)*(fu.*dfu));
    %take end point of global and local errors to find scaling for h
    lex(k)=le(end);
    gex(k)=ge(end);
    k=k+1;
    %plot for h=0.01 & h=0.02
    if h==0.01||h==0.02
        plot(0:h:0.4,le)
        title(['Local error in Euler method, h=',num2str(h)])
        xlabel('t')
        ylabel('local error')
        print('-dpng',['euler_le_h_',num2str(h),'.png'])
        clf
        plot(0:h:0.4,ge)
        title(['Global error in Euler method, h=',num2str(h)])
        xlabel('t')
        ylabel('global error')
        print('-dpng',['euler_ge_h_',num2str(h),'.png'])
    end
    clf
    clear x ge fu dfu i le
end
%plot local/global errors against h
plot([0.0025,0.005,0.01,0.02,0.04,0.08],lex,'r')
title('local error for Euler Method for changing h')
xlabel('h')
ylabel('local error')
print('-dpng',['euler_le_scaling','.png']);
clf
plot([0.0025,0.005,0.01,0.02,0.04,0.08],gex,'r')
title('global error for Euler Method for changing h')
xlabel('h')
ylabel('global error')
print('-dpng',['euler_ge_scaling','.png']);
clf
loglog([0.0025,0.005,0.01,0.02,0.04,0.08],gex,'r')
title('ln global error for Euler Method for changing h')
xlabel('ln h')
ylabel('ln global error')
print('-dpng',['euler_ln_ge_scaling','.png']);