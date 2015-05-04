function[zout]=RKq6_dzdt(t,x,y,z)
    %I know t(i) argument is not needed in the following expressions, as
    %dx/dt and dy/dt are not defined with t, but it means the code can
    %be used for other functions of t, x and y with minimal changes.
    a=0.2;
    b=0.2;
    c=5.7;

zout=b+z*(x-c);
end