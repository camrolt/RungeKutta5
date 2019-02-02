clear all
kmax = 100000;

h = .01;
y0 = 1;

y = zeros(1, kmax);
y(1) = y0;
s = 0:h:h*(kmax-1);
 
%y(k+1) = RungeKutta5(@testfunction, s, y0);
y = RungeKutta5(@testfunction, s, y0);

%y = feval(@testfunction, s);
[sode, yode] = ode45(@testfunction, s, y0);

figure
hold on
plot(s, y) %Numerical Solution
plot(sode, yode) %Numerical Solution
plot(s, s.^2 + 1) %Theoretical Solution