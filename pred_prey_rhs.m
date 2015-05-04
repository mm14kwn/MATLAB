function z_out=pred_prey_rhs(z)
p=zeros(1,2);
p(1)=-(z(1))*(2-z(2));
p(2)=(z(2))*(1-2*(z(1)));
z_out=p;