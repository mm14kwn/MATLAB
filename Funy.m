function[Fy]=Funy(t,k,i,j,add,pos)
Fy=(k(j)*(positioni(t,1,add,i,pos)-positionj(t,1,add,j,pos))/(((positioni(t,1,add,i,pos)-positionj(t,1,add,j,pos))^2)+((positioni(t,0,add,i,pos)-positionj(t,0,add,j,pos))^2)));
end