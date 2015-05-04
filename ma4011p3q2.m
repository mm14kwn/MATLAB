N=2.^(2:10); %define N
for j=1:9
   err=ma4011p3q1(N(j),0,1)-ma4011p3_exact(N(j),0,1); %calculate vector of errors at all xi for this value of N
   E(j)=max(abs(err)); % define E as maximum of modulus of this error
   clear err
end
loglog(N,E) %plot N against E
hold on
loglog(N,N.^-1,'c') %plot different values of N^-a for scale
loglog(N,N.^-2,'r')
loglog(N,N.^-3,'g')
title({'Log-Log plot of N against E(N) (in blue)';'For scale: Cyan N^{-1}, Red N^{-2}, Green N^{-3}'})
xlabel('N')
ylabel('E(N)')
hold off