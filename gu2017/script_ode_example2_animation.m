%% ODE Example #2 with animation
% T. Fitzgerald
%
% This example integrates the equations of a pendulum on a cart.
%
% $$ \left[\begin{array}{cc} M+m & -ml\cos\theta\\-m l \cos\theta & ml^2 \end{array} \right]\left[\begin{array}{c} \ddot x \\ \ddot \theta \end{array} \right]= \left[\begin{array}{c} - m l \dot \theta^2 \sin\theta + f(t) \\m g l \sin\theta\end{array} \right] $$
%
% It follows the same form as the previous example, but also includes an
% animation of the simulation.

clear all; close all; clc

%% Save movie to an external file?
% Decide to save the animation to an external movie or not
flag_save = 0;

%% Define the parameters of the system
m = 5;    %[kg]
M = 15;   %[kg]
g = 9.81; %[m/s^2]
l = 0.5;  %[m]

%%
% Define the forcing function
f = @(t,z) 0;

%%
% Define the initial conditions, where the state-vector is
% |z = [ x; theta; x dot; theta dot]|
z0 = [0; -90*pi/180; 0; 0];

%%
% Define the time to integrate
tf = 10;

%% Numerically solve
% I'll get out the solver struct, instead of the point directly.  This
% sometimes gives useful info about the solution process.  Also, I'll make
% an anonymous function to effectively pass all the parameters into the
% ode function.
sol = ode45( @(t,z) ode_pendulumcart(t, z, M, m, g, l, f),...
    [0,tf], z0);

%% Animate the solution
% I'll initiate the blank figure window, and set the aspect ratio to 1:1
% (that way that circles will be round, etc.)
fig = figure();

% make the figure's background white
set(fig, 'Color', 'w');

% make the figure big on the screen
set(fig, 'units', 'normalized', 'outerposition', [0 0 0.8 0.8])

% setup the axes
ax = axes('Parent', fig);
hold(ax, 'on');
set(ax,'DataAspectRatio', [1,1,1]);
xlim([-1,1]*l*1.5);
ylim([-1.1, 0.3]*l);
xlabel('x [m]');
ylabel('y [m]');


%%
% Setup functions to take the states and figure out locations of the
% various bodies in the figure.  I'll start with the rectangle for the
% cart.
cart_width  = l/3;
cart_height = l/6;

cart_x_pts = @(z) [ -1; +1; +1; -1]*cart_width/2 + z(1);
cart_y_pts = [ -1; -1; +1; +1]*cart_height/2;

h_cart = patch( cart_x_pts(z0), cart_y_pts, 'r', 'Parent', ax);

%%
% Next I'll add a line for the pendulum
rod_x_pts = @(z) [ z(1); z(1)-l*sin(z(2))];
rod_y_pts = @(z) [    0; l*cos(z(2))];

h_rod = line(ax, ...
    'XData', rod_x_pts(z0), ...
    'YData', rod_y_pts(z0),...
    'color', 'k', 'LineWidth', 4);

%%
% Now, I'll add the lumped mass of the pendulum.  This is the same as the
% second point of the pendulum's rod, but I need to make a second graphics
% object so I can color/size it differently.  I could use a whole point of
% points to make a circle polygon.  Instead, I'll just use a single marker
% point, make it big, and fill it in with color.

pendulum_x_pt = @(z) z(1)-l*sin(z(2));
pendulum_y_pt = @(z) l*cos(z(2));

h_pend = plot( pendulum_x_pt(z0), pendulum_y_pt(z0), ...
    'LineStyle', 'none', ...
    'marker', 'o', ...
    'MarkerSize', 20,...
    'MarkerEdgeColor', [0.1, 0.7, 0.3], ...
    'MarkerFaceColor', [0.1, 0.7, 0.3] ) ;

%%
% I'll set the title to have time info.
h_title = title(ax, 'Time = 0.0 [s]', 'FontWeight', 'bold', 'FontSize', 12);

%%
% To make the animation, I'll step through the time results, and update the
% various graphics objects.  I'll add a very short pause at the end of the
% loop to ensure that the graphics buffer is flushed.

nframes = 100;
time = linspace(0,tf,nframes);

%% 
% initialize the output file
if flag_save
    V = VideoWriter('cart_sim', 'MPEG-4');
    fps = nframes/tf;
    V.FrameRate = fps;
    open(V)
end

for i = 1:nframes
    
    z = deval( sol, time(i) );
    
    % update the title
    set( h_title, 'String', sprintf('Time = %4.1f [s]', time(i)) );
    
    % update the cart
    set( h_cart, 'XData', cart_x_pts(z) );
    
    % update the rod
    set( h_rod, 'XData', rod_x_pts(z), ...
        'YData', rod_y_pts(z) );
    
    % update the point-mass location
    set(h_pend, 'XData', pendulum_x_pt(z), ...
        'YData', pendulum_y_pt(z));
    
    % pause so we can actually watch the animation
    pause(0.01)
    
    % caputure an image of the figure, and write it to the movie
    if flag_save
        drawnow
        writeVideo( V, getframe(fig) );
    end
    
end

% close the movie file
if flag_save
    close(V)
end
