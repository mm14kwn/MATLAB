function [x]=ma3012p2q1(A,b)
n=length(b)
for i=1:n
    l(i)=i
    s(i)=max(abs(A(i,:)))
end
for k=1:n-1
    for v=k:n
        for u=1:n
            w(u)=abs(A(u,k))/s(u)
            y(u)=abs(A(v,k))/s(v)
            if w(u)>=y(u)
                p=u
                break
                break
            end
        end
    end
    if p~=k
        h=l(p)
        l(p)=l(k)
        l(k)=l(p)
    end
    if A(l(k),k)==0
        error('Matrix A is singular')
    end
    for i=k+1:n
        m(l(i),k)=A(l(i),k)/A(l(k),k);  
        for j=k+1:n
            A(l(i),j)=A(l(i),j)-(m(l(i),k)*A(l(k),j));
        end
        b(l(i))=b(l(i))-(m(l(i),k)*b(l(k)));
    end
end
x=kwn_bck_sub(A,b);
end
    