function z=dotproduct(x,y)
a=length(x);
b=length(y);
if a~=b
    error('Please enter two vectors of the same size')
    return
end
z=0;                                                           
for i=1:a;                                                    
z=z+x(i)*y(i);
end
end