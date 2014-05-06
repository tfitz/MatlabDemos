%% animate the pendulum
% The pendulum is subjected to a small initial disturbance.

%% clear the workspace
clear all;
close all;
clc


%% Define the parameters of the system
g =   9.81 ; % [m/s^2]
L =   1.   ; % [m]
k = 400.   ; % [N/m]
m =   2.   ; % [kg]
a =   L/4  ; % [m]

tf = 10; % [s]
nframes = 60;

%%
% Animate the solution?
flag.animate = false;
flag.output.gif = false;

%% define the nonlinear equation of motion
% Here "x" is a stand-in for the two states of the system $\theta$ and 
% $\dot{theta}$.
eom = @(t,x,f) [...
            x(2,:); ...
            (f(t,x) - (a^2*k*cos(x(1,:)) - m*g*L).*sin(x(1,:)))/(m*L^2)...
            ];
f = @(t,x) 0;


%% Define the initial conditions
% Since we are considering an impulse response, the initial conditions are
% set as a scaled initial velocity.
x0 = [0;
      1/(m*L^2)];
  
  
%% Numerical solve the ODE
sol = ode45(@(t,x) eom(t,x,f), [0, tf], x0);
time = linspace(0, tf, 1000);
x = deval(sol,time)';

%% Plot the results
fig = figure();
ax = axes('parent', fig);
hold(ax, 'on');
plot(ax, time, x(:,1)*180/pi, 'b-')
xlabel(ax, 'time [s]', 'FontSize', 12)
ylabel(ax, 'Angular position \theta [deg]', 'FontSize', 12)


%% Animate the results
if( flag.animate )
    
    time_animation = linspace(0,tf,nframes);
    
    fig = figure();
    set(fig, 'position', [293, 337, 1139, 611]);
    set(fig, 'color', 'w');
    
    %%
    % setup subplot 1
    ax1 = subplot(1,3,1:2,'Parent',fig);
    set(ax1, 'DataAspectRatio', [1 1 1]);
    xlabel(ax1, 'x/L position', 'FontSize', 12)
    ylabel(ax1, 'y/L position', 'FontSize', 12)
    xlim(ax1, [-1.1, 1.1]);
    ylim(ax1, [-0.5, 1.1]);
    
    h_title = title(ax1,'Time = 0 [s]');
    
    h_base = rectangle('Parent', ax1,...
        'Position', [-1/4, -1/16, 1/2, 1/16]);
    
    h_line = line('Parent', ax1,...
        'XData', [0, L*sin(x0(1))]/L,...
        'YData', [0, L*cos(x0(1))]/L,...
        'LineWidth',3);
    
    %%
    % Let's assume that the spring remains fixed at one end.
    h_spring = line('Parent', ax1, ...
        'XData', [ a*sin(x0(1)), 0.5*L]/L, ...
        'YData', [ a, a]/L, ...
        'LineStyle', '--', ...
        'Color', 'b', ...
        'LineWidth', 3);
    
    h_base_spring = rectangle('Parent', ax1,...
        'Position',[0.5*L, a-a/4, L/16, a/2]/L);
    
    %%
    % Let's take the cheap way out, and draw a circle just as a marker point
    h_tip = patch('Parent', ax1,...
        'XData', sin(x0(1)),...
        'YData', cos(x0(1)),...
        'CData', 1,...
        'Marker','o',...
        'MarkerFaceColor','r',...
        'EdgeColor','r',...
        'MarkerSize',24);
    
    %%
    % Let's setup subplot 2
    ax2 = subplot(1,3,3,'Parent',fig);
    hold(ax2, 'on');
    plot(ax2, time, x(:,1)*180/pi, 'b-','LineWidth', 2)
    xlabel(ax2, 'time [s]', 'FontSize', 12)
    ylabel(ax2, 'Angular position  \theta [deg]', 'FontSize', 12)
    
    h_trace =  patch('Parent', ax2,...
        'XData', 0,...
        'YData', x0(1)*180/pi,...
        'CData', 1,...
        'Marker', 'o',...
        'MarkerFaceColor', 'm',...
        'EdgeColor', 'm',...
        'MarkerSize', 10);
    
    
    %%
    % Now let's loop over each timestep and grab the frame.  We'll save each
    % frame into an animated GIF.
       
    if( flag.output.gif )
        drawnow;
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        outfile_name = sprintf('pendulum_k=%d.gif',k);
        imwrite(imind,cm,outfile_name,'gif','DelayTime',0,'loopcount',inf);
    end
    
    for i = 1:nframes
        
        theta = deval(sol, time_animation(i), 1);

        set(h_line,...
            'XData', [0, sin(theta)],...
            'YData', [0, cos(theta)]);
        
        set(h_tip,...
            'XData', sin(theta), ...
            'YData', cos(theta) );
        
        set(h_spring,...
            'XData', [ a*sin(theta), 0.5*L]/L, ...
            'YData', [ a*cos(theta), a]/L);
        
        set(h_trace,...
            'XData', time_animation(i), ...
            'YData', theta*180/pi );
        
        set(h_title,'String',sprintf('Time = %6.2f [s]',time_animation(i)));

        if( flag.output.gif )
            drawnow;
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,outfile_name,'gif','DelayTime',0,'writemode','append');
        else
            pause(0.05);
        end

    end   
          
          
end
