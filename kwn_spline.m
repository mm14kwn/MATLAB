function s=kwn_spline(t,t_in,y,y2drv)
if size(t_in)~=size(y)
    error('please enter t,y,y2drv the same size')
end
if size(y)~=size(y2drv)
    error('please enter t,y,y2drv the same size')
end
if t<min(t_in)||t>max(t_in)
    error('please enter t within the bounds of t0:tn')
end
if t==t_in(1)
    k=1;
end
if t==t_in(length(t_in))
    k=length(t_in)-1;
end
for i=1:length(t_in)-1
    if t_in(i)>=t_in(i+1)
        error('please enter t_in as a vector of distinct, ordered values')
    end
    if t_in(i)<=t && t<t_in(i+1)
        k=i;
        break
    end
end
h=t_in(k+1)-t_in(k);
s=y(k)+((1/h)*(y(k+1)-y(k))-(h/6)*(y2drv(k+1)+2*(y2drv(k))))*(t-t_in(k))+(1/2)*(y2drv(k))*((t-t_in(k))^2)+(1/(6*h))*(y2drv(k+1)-y2drv(k))*((t-t_in(k))^3);
end


