function out=f_ellipse(u,v,a,b,c)
out=a*b*c*cos(u)*sin(v)*sin(u)*sin(v)*cos(v)*sin(v)*sqrt(((((b^2)*(c^2)*((cos(u))^2))+((a^2)*(c^2)*((sin(u))^2)))*((sin(v))^2))+((a^2)*(b^2)*((cos(v))^2)));
end