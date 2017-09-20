%% ODE Example #2
% T. Fitzgerald
%
% This example integrates the equations of a pendulum on a cart.  The
% equations are put implemented in the following form
%
% $$ \left[\begin{array}{cc} M+m & -ml\cos\theta\\-m l \cos\theta & ml^2 \end{array} \right]\left[\begin{array}{c} \ddot x \\ \ddot \theta \end{array} \right]= \left[\begin{array}{c} - m l \dot \theta^2 \sin\theta + f(t) \\m g l \sin\theta\end{array} \right] $$
%
% This gives me a way to sidestep having to put the equations into explicit
% first-order form.  I can invert the matrix on the right, and then build
% state-space form.

%%

clear all
close all
clc

%% Define the parameters of the system
m = 5;    %[kg]
M = 15;   %[kg]
g = 9.81; %[m/s^2]
l = 0.5;  %[m]

%%
% Define the forcing function
f = @(t,y) 0;

%%
% Define the initial conditions, where the state-vector is
% |y = [ x; theta; x dot; theta dot]|
y0 = [0; -90*pi/180; 0; 0];

%%
% Define the time to integrate
tf = 10;

%% Numerically solve
% I'll get out the solver struct, instead of the point directly.  This
% sometimes gives useful info about the solution process.  Also, I'll make
% an anonymous function to effectively pass all the parameters into the
% ode function.
sol = ode45( @(t,y) ode_pendulumcart(t, y, M, m, g, l, f),...
    [0,tf], y0);

%% Plot
% To generate a plot, we'll need to use some other functions.  The data
% contained in |sol| is probably very coarse, so I'll want to smooth it up
% to make a nice looking plot.

t     = linspace(0,tf,300);
x     = deval( sol, t, 1);
theta = deval( sol, t, 2);

figure
subplot(2,1,1)
plot( t, x/l, 'LineWidth', 1.5)
ylabel('Cart position x/l')
grid on

subplot(2,1,2)
plot(t, theta*180/pi, 'LineWidth', 1.5)
ylabel('Pendulum position \theta [deg]')
xlabel('Time t [s]')
grid on
