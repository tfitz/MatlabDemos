%% Mode Shape Animation
% This demonstration uses matlab to calculate the eigenvalues and vectors
% of the two-mass system.  It then uses this as initial
% conditions to simulate the free response of the system.  A brief summary
% of the problem is available at
% <http://fitzgeraldlab.wordpress.com/2013/12/19/example-of-animating-modeshapes>

clear all; close all; clc

%% Choose the graphic outputs
% Animate the solution?
flag.animate = true;
flag.output.gif = true;

%% Initial conditions
% Pick the initial conditions, as being scaled shapes from the modeshapes
%   y0 = alpha1*"mode shape 1" + alpha2*"mode shape 2"
alpha1 = 0;
alpha2 = 0.5;

%% Define system parameters
% in SI units
m1 = 1;
m2 = 1;
k1 = 10;
k2 = 15;
L1 = 1;
L2 = 1;
c1 = 0;

%%
% The solution time
tf = 4;

%% Build the system matrices
M = diag([m1 m2])

a = 1+L1/L2;
b = -L1/L2;
K = [ k1*a^2, k1*a*b;
    k1*a*b, k2+k1*b^2]

C = [c1, 0;
    0, 0]

%% Solve the algebraic eigenvalue problem
% This solves the undamped eigenvalue problem
[Phi,lambda] = eig(K,M);

%%
% The $\omega$ values are
omega = sqrt(diag(lambda))

%%
% Let's make the columns of $[\Phi]$ look like the common way it does in
% class:
for i = 1:2
    Phi(:,i) = Phi(:,i)/Phi(1,i);
end
Phi

%% The initial conditions
% Use the eigenvectors to build the initial conditions
y0 = alpha1*Phi(:,1) + alpha2*Phi(:,2);

%% Numerical simulation
% Setup the state-space form of the equations
A = [ zeros(2,2), eye(2);
    -M\K    , -M\C ];
z0 = [y0; zeros(2,1)];

sol = ode45(@(t,z) A*z, [0 tf], z0);
t = linspace(0,tf,1000);
z = deval(sol,t);

%% Plot the response as a function of time
figure
plot(t,z(1,:),'b','Displayname','y_1');
hold on
plot(t,z(2,:),'g--','Displayname','y_2');
legend('show');
xlabel('time t');
ylabel('Position of each mass');

%% Animate the solution
if( flag.animate )

    nframes = 300;
    time = linspace(0,tf,nframes);
    z = deval(sol,time,[1,2]);
    Yref = L1+L2;

    fig = figure('color','w');
    ax = axes('parent',fig);
    hold(ax,'on');
    y3 = @(y1,y2) a*y1 + b*y2;
    ymin = -2;
    ymax = 2;
    ylim(ax,[ymin ymax]);
    xlim(ax,[-0.1 1.1]*(L1+L2));
    set(ax,'XTick',[0,L1,L1+L2])
    set(ax,'XTickLabel',{'0','L1','L1+L2'});
    h_title = title(ax,'time = 0');
    ylabel(ax,'Normalized position y_i/(L_1+L_2)');

    %%
    % init. the static line
    h_static = plot(ax,[0, L1+L2],[0 0]/Yref,'--','LineWidth',3,'Color',0.8*[1 1 1]);

    %%
    % init. the moving line
    h_bar = plot(ax,[0,L1+L2],[y3(y0(1),y0(2)),y0(2)]/Yref,'-','LineWidth',2,'Color','k');

    %%
    % init the springs damper
    if( c1 > 0 )
        h_c1 = plot(ax, [L1,L1],[ymin,y0(1)]/Yref,':','LineWidth',3,'Color','r');
    end
    h_k1 = plot(ax, [0,0],[ymin,y3(y0(1),y0(2))]/Yref,':','LineWidth',3,'Color',[205, 78, 212]/255);
    h_k2 = plot(ax, [L1+L2,L1+L2],[ymin,y0(2)]/Yref,':','LineWidth',3,'Color',[205, 78, 212]/255);

    %%
    % init. the dots
    h_y1 = plot(ax,L1,y0(1)/Yref,'Marker','o','color','b','MarkerFaceColor','b','MarkerSize',14);
    h_y2 = plot(ax,L1+L2,y0(2)/Yref,'Marker','o','color','g','MarkerFaceColor','g','MarkerSize',14);

    if( flag.output.gif )
        drawnow;
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        outfile_name = 'TwoParticleTimeHistory.gif';
        imwrite(imind,cm,outfile_name,'gif','DelayTime',0,'loopcount',inf);
    end

    for i = 1:nframes

        Y1 = z(1,i);
        Y2 = z(2,i);
        Y3 = y3(Y1,Y2);

        set(h_bar,'YData',[Y3, Y2]/Yref);
        set(h_y1,'YData',Y1/Yref);
        set(h_y2,'YData',Y2/Yref);

        if( c1 > 0 )
            set(h_c1,'YData',[ymin,Y1]);
        end
        set(h_k1,'YData',[ymin,Y3/Yref]);
        set(h_k2,'YData',[ymin,Y2/Yref]);

        set(h_title,'String',sprintf('time = %.3f',time(i)));

        if( flag.output.gif )
            drawnow;
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,outfile_name,'gif','DelayTime',0,'writemode','append');
        else
            pause(0.04);
        end

    end

end
