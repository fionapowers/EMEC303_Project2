%% HW 4
%% Heun's Method: Problem 1
clear; clc

%% Inputs
%h = height of water above the outlet
%Do = diameter of the circular outlet
%Dt = diameter of the cylindrical tank
%g = gravity m/s*^2
%Assume Dt>>Do
%Assume P (Pressure) = constant
%Bernoulli's equation:
    % ((vt)*^2/2) + g*(zt) + (Patm/rho) = ((vo)*^2/2) + g*zo + (Patm/rho)
    
%h = zt-zo
%vt = 0 m/2

%vo = sqrt(2*g*h)
%-At*(dh/dt) = Ao*vo

%(dh/dt) = -(Do*^2/Dt*^2)*sqrt(2*g*h)

%% Inputs
h = 3; 
Dt = 1; 
Do = 0.01;
g = 9.8;
yo = 1;
xl = 0;
xh = 100;

% Define Function
%f = @(x,y) 2*cos(x)+0.2; 
f = @(x,y) -(Do.^2/Dt.^2)*sqrt(2*g*y);

% Initialize Arrays
x = 0;%[xl:h:xh];
y = zeros(1,10000); % heuns

% Set Initial Conditions
y(1) = yo;
yold = yo;
i=1;
% Iterate Over Space
while yold > 0 %i = 1:10000%length(x)-1
    yold = y(i);
    %Heun
    ystar = y(i) + f(x,y(i))*h;
    y(i+1) = y(i) + (f(x,y(i)) + f(x+h,ystar))*h/2; %correction
    x = x + h;
    i = i+1;
end

% Plot and compare 
figure(1); clf(1)
plot(0:h:x,y(1:i),'-*')
hold on
xlabel('Time (s)')
ylabel('Height (m)')
title('Plot 1: Tank Problem')


%% Problem 2 RK-4 Method & Analytical Method
% EMEC303 RK4 vs Analytic

%Inputs
h = 0.0015;
Lx = .05;
N = round(Lx/h); % numebe of steps in loop
yo = 0; % initial condition


%Define Function
f = @(x,y) -1000*y+3000-2000.^(-x);
fA = @(x) 3-0.998.^(-1000*x)-2.002.^(-x);

% Est Arrays/Preallocate
x = zeros(1,N);
y = zeros(1,N);

% Set Initial Condition
y(1) = yo; % set Initial condition after preallocate

% Numerical Approx-RK4

for i=1:N
    
    % New x value
    x(i+1)=x(i)+h;

    % New y value using 4th order Runge-Kutta Method
    k1=f(x(i)      ,y(i)         )';
    k2=f(x(i)+0.5*h,y(i)+0.5*k1*h)';
    k3=f(x(i)+0.5*h,y(i)+0.5*k2*h)';
    k4=f(x(i)+    h,y(i)+    k3*h)';
    y(i+1)=y(i)+1/6*(k1+2*k2+2*k3+k4)*h;
    
end

% Plot
figure(1); clf(1)
plot(x,y)
hold on

plot(x, fA(x), '*')
legend('RK-4', 'Analytic')
hold on

    

%% Problem 2 Euler's Method 

% Euler's Method
% Inputs

h = 0.000015;
yo = 1; %initial Condition
tl = 0;
th = 1;

% Establish domain
t = [tl:h:th]; %low and high step points, t(i+1) = t(i)+h;
y = zeros(1, length(t)); %preallocate

% Set Initial Condition

y(1) = yo; %remember no 0th interval in MATLAB/ after preallocate

% Define Function
f = @(x,y) -1000*y+3000-2000.^(-x);

% Iterate (loop over space/time)
for i = 1:length(t)-1
    % t(i+1) = t(i)+h; %intependent variable 
    y(i+1) = y(i)+f(y(i),t(i)).*h; %finish eq. %dependent variable
end

% Plot Solution
figure(1); clf(1);
plot(t,y,'-*')
hold on
xlabel('time')
ylabel('y')
title('Euler''s Method')