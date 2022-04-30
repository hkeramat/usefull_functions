x = pi*[0:.5:2]; 
y = [0  1  0 -1  0  1  0; 
     1  0  1  0 -1  0  1];
pp = spline(x,y);
yy = ppval(pp, linspace(0,2*pi,101));
plot(yy(1,:),yy(2,:),'-b',y(1,2:5),y(2,2:5),'or')
axis equal
