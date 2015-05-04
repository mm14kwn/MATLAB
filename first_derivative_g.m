clear all;                      %       clears all the variables from memory
n=1000;                           % number of points used in the subdivision of the interval [0,1]
h =1/n;                         % set up the difference step size
a=(2*pi)/3;
x = linspace(0,a,n+1);          %  create the vector x = (0; h; 2h;… ; 1)
x_plus_h = x + h.*ones(1,n+1); %   create the vector x plus h= (h; 2h; 3h;… ; 1; 1 + h)
x_minus_h = x - h.*ones(1,n+1); % create the vector x minus h= (-h; 0; h; 2h;…; 1- h)
g_prime = (g(x_plus_h)-g(x_minus_h))./(2*h);    % get the divided difference 
plot(x,g_prime)                 %      plot the derivative