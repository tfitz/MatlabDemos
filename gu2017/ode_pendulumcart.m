function dy = ode_pendulumcart( t, y, M, m, g, l, f)

%% Equations of motion
%
% $$ (M+m)\ddot x - m l \ddot \theta \cos \theta + m l \dot{\theta}^2 \sin\theta = f(t) $$
%
% $$ - m l \cos\theta \ddot x  + m l^2 \ddot \theta - m g l \sin\theta = 0 $$
%
% Instead of explicitly solving this into first order form, I'll just
% numerically solve it.

%%
% Unpack the state vector
x     = y(1);
theta = y(2);
xd    = y(3);
thetad= y(4);

%%
% I'll group the equations into the following form
%
% $$ \left[\begin{array}{cc} M+m & -ml\cos\theta\\-m l \cos\theta & ml^2 \end{array} \right]\left[\begin{array}{c} \ddot x \\ \ddot \theta \end{array} \right]= \left[\begin{array}{c} - m l \dot \theta^2 \sin\theta + f(t) \\m g l \sin\theta\end{array} \right] $$
%

H = [           M+m, -m*l*cos(theta);
    -m*l*cos(theta),          m*l^2];

r = [-m*l*thetad^2*sin(theta) + f(t,y);
    +m*l*g*sin(theta)];

%%
% Pack up the output $\dot y$
dy = [xd;
    thetad;
    H\r];


