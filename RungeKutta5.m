function [F] = RungeKutta5(funct, tspan, f0) %call RungeKutta5(@func, double (like 1:20 sec), double[x, y, ...])
% RUNGEKUTTA5 5th order Runge-Kutta, works with scalers and vectors
%input -> h is integration small inc
%input -> s is <x,y> between points (can be vec)
%input -> func is function name where y' = f(y)
%output <- y1 = y^(n+1) next step
F = zeros(1, length(tspan));
F(1) = f0;

for i = 1:length(tspan)-1
h = tspan(i+1) - tspan(i); %integration const

f = funct(tspan(i), 0);
k1 = h*f;
k2 = h*funct(tspan(i), F(i) + (k1/3) + (h*f*k1/18)); 
k3 = h*funct(tspan(i), F(i) - (152*k1/125) + (252*k2/125) - (44*h*f*k1/125));
k4 = h*funct(tspan(i), F(i) + (19*k1/2) - (72*k2/7) + (25*k3/14) +(5*h*f*k1/2));

F(i+1) = F(i) + (5*k1/48) + (27*k2/56) + (125*k3/336) + (k4/24);
end


end

