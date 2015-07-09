function plot1(p,e,t,ustore,vstore,hstore,ti)
  pdeplot(p,e,t,'zdata',hstore(:,ti),'xydata',[ustore(:,ti),vstore(:,ti)],'contour','on')
hold on
  axis([min(p(1,:)) max(p(1,:)) min(p(2,:)) max(p(2,:)) -1 1])
  title(int2str(ti))
hold off
end
