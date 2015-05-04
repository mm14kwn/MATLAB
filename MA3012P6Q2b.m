clear all
k=1;
for h=[0.0025,0.005,0.01,0.02,0.04,0.08,0.16]
    %set initial conditions
    x(1)=0.1;
    ge(1)=0;
    fu(1)=10/9;
    dfu(1)=100/81;
    d2fu(1)=2000/729;
    d3fu(1)=20000/2187;
    d4fu(1)=800000/19683;
    i=2;
    %calculate values at each point
    for t=h:h:0.4
        xhold=x(i-1)+h*(1/(1-x(i-1)));
        x(i)=x(i-1)+(h/2)*((1/(1-x(i-1)))+(1/(1-xhold)));
        fu(i)=1/(1-x(i));
        dfu(i)=1/((x(i)-1)^2);
        d2fu(i)=2/((x(i)-1)^3);
        d3fu(i)=6/((x(i)-1)^4);
        d4fu(i)=24/((x(i)-1)^5);
        ge(i)=abs((1-sqrt(0.81-2*t))-x(i));
        i=i+1;
    end
    %calculate local error
    le=abs((h^2)*(1/2)*(fu.*dfu))+((1/4)*(h^3)*(fu.^2).*d2fu)+((1/12)*(h^4)*(fu.^3).*d3fu)+((1/48)*(h^5)*(fu.^4).*d4fu);
    %take middle point of global and local errors to find scaling for h
    lex(k)=le(ceil(end/2));
    gex(k)=ge(ceil(end/2));
    k=k+1;
    %plot for h=0.01 & h=0.02
    if h==0.01||h==0.02
        plot(0:h:0.4,le);
        title(['Local error in Improved Euler method, h=',num2str(h)]);
        xlabel('t');
        ylabel('local error');
        print('-dpng',['imp_euler_le_h_',num2str(h),'.png']);
        clf
        plot(0:h:0.4,ge);
        title(['Global error in Improved Euler method, h=',num2str(h)]);
        xlabel('t');
        ylabel('global error');
        print('-dpng',['imp_euler_ge_h_',num2str(h),'.png']);
    end
    clf
    clear x ge fu dfu d2fu d3fu d4fu i le
end
%plot local/global errors against h, and loglog to determine scaling
plot([0.0025,0.005,0.01,0.02,0.04,0.08,0.16],lex,'r')
title('Midpoint local error for Improved Euler Method for changing h')
xlabel('h')
ylabel('Midpoint local error')
print('-dpng',['imp_euler_le_scaling','.png']);
clf
loglog([0.0025,0.005,0.01,0.02,0.04,0.08,0.16],lex,'r')
title('Ln Midpoint local error for Improved Euler Method for changing h')
xlabel('ln h')
ylabel('ln Midpoint local error')
print('-dpng',['imp_euler_ln_le_scaling','.png']);
clf
plot([0.0025,0.005,0.01,0.02,0.04,0.08,0.16],gex,'r')
title('Midpoint global error for Improved Euler Method for changing h')
xlabel('h')
ylabel('Midpoint global error')
print('-dpng',['imp_euler_ge_scaling','.png']);
clf
loglog([0.0025,0.005,0.01,0.02,0.04,0.08,0.16],gex,'r')
title('ln Midpoint global error for Improved Euler Method for changing h')
xlabel('ln h')
ylabel('ln Midpoint global error')
print('-dpng',['imp_euler_ln_ge_scaling','.png']);