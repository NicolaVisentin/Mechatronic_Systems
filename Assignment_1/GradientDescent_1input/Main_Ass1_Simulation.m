% Nicola Visentin, 20/10/2024, Politecnico di Milano

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       ASSIGNMENT 1-2: SHUTTLE                         %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                            Simulation                                 %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Animation settings

FPS=60;
duration=5;

% set(0, 'DefaultAxesFontSize', 14);   
% set(0, 'DefaultLineLineWidth', 2);    
% set(0, 'DefaultAxesLineWidth', 1.5); 
% set(0, 'DefaultAxesFontSize', 14);  
% set(0, 'DefaultLineLineWidth', 2);  
% set(0, 'DefaultAxesLineWidth', 1.5); 

set(0, 'DefaultAxesFontSize', 'default');   
set(0, 'DefaultLineLineWidth', 'default');    
set(0, 'DefaultAxesLineWidth', 'default'); 
set(0, 'DefaultAxesFontSize', 'default');  
set(0, 'DefaultLineLineWidth', 'default');  
set(0, 'DefaultAxesLineWidth', 'default'); 

%% Data and parameters definition

% shuttle data

Cd=0.5;             % drag coefficient [-]
Cl=0.2;             % lift coefficient [-]
alpha=0.00005;      % fuel consumption constant [(kg/s)/N]
g=9.81;             % gravity [m/s^2]
S=2;                % areodynamic surface (hyp.: S=1)
rho=1;              % air density (hyp.: rho=1=const)

param.Cd=Cd;
param.Cl=Cl;
param.alpha=alpha;
param.g=g;
param.S=S;
param.rho=rho;

% simulation data

t0=0;           % initial time [s]
tf=110;         % final time   [s]

h_i=100;        % initial altitude [m]
v_i=100;        % initial velocity [m/s]
m_i=20000;      % initial mass (shuttle + fuel) [kg]
gamma_i=pi/3;   % initial angle [rad]

x_in=[h_i v_i m_i gamma_i]';

% control definition

% N=10001;           
% t_u=linspace(t0,tf,N);
% 
% evolution=cos(2*pi/tf*t_u).*exp(-0.01*t_u)+0.1*t_u;
% u=m_i*g*(ones(1,N)+evolution);

load('optimal_u.mat','u')
N=length(u);
t_u=linspace(t0,tf,N);

%% Simulation

[t_x,x]=ode45(@(t,x) EquationOfMotion(t,x,u,t_u,param),[t0 tf],x_in);

%% Animation

% setup

N_anim=FPS*duration;              % number of frames ( = number of data points)
t_anim=linspace(t0,tf,N_anim);
x_anim=interp1(t_x,x,t_anim);
u_anim=interp1(t_u,u,t_anim);

h_anim=x_anim(:,1);                       % height
v_anim=x_anim(:,2);                       % speed
m_anim=x_anim(:,3);
gamma_anim=x_anim(:,4);                   % inclination

len_sh=(0.1*(max(h_anim)-min(h_anim)));       % "shuttle" (segment) length
len_fl=u_anim*len_sh/max(u_anim);             % "flame" length: visualisation of u

yaxis_sh=[min(h_anim)-2*len_sh, max(h_anim)+2*len_sh];  % fix ylim
xaxis_sh=yaxis_sh-mean(yaxis_sh);                       % dynamic xlim
yaxis_tr=[min(h_anim)-2*len_sh, max(h_anim)+2*len_sh];  % fix ylim
xaxis_tr=yaxis_sh-mean(yaxis_sh);                       % dynamic xlim

% open figure

AnimFig=figure;
AnimFig.WindowState='maximized';

subplot(2,4,3)   % shuttle
hold on
pl_sh=plot(nan,nan,'k','LineWidth',2);
pl_fl=plot(nan,nan,'r','LineWidth',4);
pl_shtr=plot(nan,nan,'k--','LineWidth',1);
grid on
box on
ax_sh=gca;
ax_sh.XGrid='off';
ax_sh.XTickLabel=[];
ylabel('h [m]','Interpreter','latex')
title('Animation')
xlim(xaxis_sh)
ylim(yaxis_sh)
pl_ground=fill([xaxis_sh(1) xaxis_sh(1) xaxis_sh(2) xaxis_sh(2)],[0 yaxis_sh(1) yaxis_sh(1) 0],[0.8 0.8 0.8],'FaceAlpha',0.9);

subplot(2,4,[1 2])  % trajectory
hold on
pl_tr=plot(nan,nan,'k--','LineWidth',1);
pl_tr_mark=plot(nan,nan,'ko','MarkerFaceColor','k');
grid on
box on
ax_tr=gca;
xlabel('x [m]','Interpreter','latex')
ylabel('h [m]','Interpreter','latex')
title('Trajectory')
xlim(xaxis_tr)
ylim(yaxis_tr)
pl_ground2=fill([xaxis_tr(1) xaxis_tr(1) xaxis_tr(2) xaxis_tr(2)],[0 yaxis_tr(1) yaxis_tr(1) 0],[0.8 0.8 0.8],'FaceAlpha',0.9);

subplot(2,4,5)     % altitude
hold on
pl_h=plot(nan,nan,'b');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('h [m]','Interpreter','latex')
title('Altitude')
xlim([t0 tf]);
ylim([min(h_anim) max(h_anim)])

subplot(2,4,6)     % velocity
hold on
pl_v=plot(nan,nan,'g');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Velocity')
xlim([t0 tf])
ylim([min(v_anim) max(v_anim)])

subplot(2,4,7)     % mass
hold on
pl_m=plot(nan,nan,'c');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('m [kg]','Interpreter','latex')
title('Mass')
xlim([t0 tf])
ylim([0 max(m_anim)])

subplot(2,4,8)     % gamma
hold on
pl_gamma=plot(nan,nan,'m');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\gamma$ [deg]','Interpreter','latex')
title('Inlcination')
xlim([t0 tf])
ylim([-360 360])

subplot(2,4,4)     % control u
hold on
pl_u=plot(nan,nan,'r');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('u [N]','Interpreter','latex')
title('Thrust')
xlim([t0 tf])
ylim([min(u_anim) max(u_anim)])

% starting position

xb=0;           % back end x coordinate
yb=h_anim(1);   % back end y coordinate

xt=xb+len_sh*cos(gamma_anim(1));  % front end x coordinate
yt=yb+len_sh*sin(gamma_anim(1));  % front end y coordinate

xu=xb-len_fl(1)*cos(gamma_anim(1));  % control visulization: x "flame" coord
yu=yb-len_fl(1)*sin(gamma_anim(1));  % control visulization: y "flame" coord

set(pl_sh,'XData',[xb xt],'YData',[yb yt])
set(pl_tr_mark,'XData',xb,'YData',yb)

pause(2)

% animation

x_tr=zeros(N_anim);
x_tr(1)=xb;
for ii=2:N_anim

    % update back end coordinates

    dt_ii=t_anim(ii)-t_anim(ii-1);
    dx_ii=v_anim(ii-1)*dt_ii*cos(gamma_anim(ii-1));

    xb=xb+dx_ii;
    yb=h_anim(ii);
    x_tr(ii)=xb;
    
    % update front end coordinates

    xt=xb+len_sh*cos(gamma_anim(ii));
    yt=yb+len_sh*sin(gamma_anim(ii));

    % update flame coordinates

    xu=xb-len_fl(ii)*cos(gamma_anim(ii));
    yu=yb-len_fl(ii)*sin(gamma_anim(ii));

    % update plots
    
    xaxis_sh=xaxis_sh+dx_ii;
    xaxis_tr=[min(xaxis_sh(1),xaxis_tr(1)), max(xaxis_sh(2),xaxis_tr(2))];

    subplot(2,4,3)
    xlim(xaxis_sh)
    set(pl_sh,'XData',[xb, xt],'YData',[yb, yt])
    set(pl_fl,'XData',[xb, xu],'YData',[yb, yu])
    set(pl_ground,'XData',[xaxis_sh(1) xaxis_sh(1) xaxis_sh(2) xaxis_sh(2)],'YData',[0 yaxis_sh(1) yaxis_sh(1) 0])
    set(pl_shtr,'XData',x_tr(1:ii),'YData',h_anim(1:ii))

    subplot(2,4,[1 2])
    set(pl_tr,'XData',x_tr(1:ii),'YData',h_anim(1:ii))
    set(pl_tr_mark,'XData',x_tr(ii),'YData',h_anim(ii))
    set(pl_ground2,'XData',[xaxis_tr(1) xaxis_tr(1) xaxis_tr(2) xaxis_tr(2)],'YData',[0 yaxis_tr(1) yaxis_tr(1) 0])
    xlim(xaxis_tr)

    set(pl_h,'XData',t_anim(1:ii-1),'YData',h_anim(1:ii-1))
    set(pl_v,'XData',t_anim(1:ii-1),'YData',v_anim(1:ii-1))
    set(pl_m,'XData',t_anim(1:ii-1),'YData',m_anim(1:ii-1))
    set(pl_gamma,'XData',t_anim(1:ii-1),'YData',rad2deg(gamma_anim(1:ii-1)))
    set(pl_u,'XData',t_anim(1:ii-1),'YData',u_anim(1:ii-1))

    % pause

    drawnow
    pause(1/FPS);

    % % save GIF
    % frame = getframe(AnimFig);
    % im = frame2im(frame);
    % [A, map] = rgb2ind(im, 256);
    % if ii == 2
    %     imwrite(A, map, 'animationSimulation.gif', 'gif', 'LoopCount', inf, 'DelayTime', 0.001);
    % else
    %     imwrite(A, map, 'animationSimulation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.001);
    % end

end
% imwrite(A, map, 'animationSimulation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 2);
pause(2)

%% Plots

% plot control in time

figure;
plot(t_u,u,'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('u [N]','Interpreter','latex')
title('Control input (thrust)')
axis tight

% plot state in time

figure;

subplot(2,2,1)
plot(t_x,x(:,1),'b')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('h [m]','Interpreter','latex')
title('Altitude profile')

subplot(2,2,2)
plot(t_x,x(:,2),'g')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Velocity profile')

subplot(2,2,3)
plot(t_x,x(:,3),'c')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('m [kg]','Interpreter','latex')
title('Mass evolution')

subplot(2,2,4)
plot(t_x,rad2deg(x(:,4)),'m')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\gamma$ [deg]','Interpreter','latex')
title('Angle evolution')