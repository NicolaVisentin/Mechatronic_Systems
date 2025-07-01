clear
close all
clc

%% Animation settings

FPS=60;
duration=2;

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

%% Data and parameters

% For the system

m=0.5;          % mass of the piston [kg]
c=1;            % damping [N s/m]
k1=1;           % stiffness coefficient [N/m]
k3=0.1;         % stiffness coefficient [N/m^3]
alpha=50;       % force constant [N/A]
L0=0.01;        % inductance at reference temperature [H]
beta1=-7e-5;    % inductance coefficient [H/°C]
beta2=2e-7;     % inductance coefficient [H/°C^2]
RES=10;         % resistance [Ohm]
Ct=500;         % thermal capacity [J/°C]
h=10;           % thermal dissipation coefficient [W/°C]
Ta=20;          % ambient temperature [°C]

param.m=m;
param.c=c;
param.k1=k1;
param.k3=k3;
param.alpha=alpha;
param.L0=L0;
param.beta1=beta1;
param.beta2=beta2;
param.RES=RES;
param.Ct=Ct;
param.h=h;
param.Ta=Ta;

% For the simulation

t0=0;        % initial time [s]
tf=5;       % final time   [s]

x1_i=0;      % initial position [m]
x2_i=0;      % inital speed [m/s]
x3_i=0;      % initial current [A]
x4_i=20;     % initial temperature [°C]

x_in=[x1_i x2_i x3_i x4_i]';

% Control definition

N=100001;
t_u=linspace(t0,tf,N);

% u=0.1*sin(t_u);             % [V]
load("actual_control.mat")
u=interp1(t_sim,u,t_u);


%% Simulation

[t_x,x]=ode78(@(t,x) EquationOfMotion(t,x,u,t_u,param),[t0 tf],x_in);

%% Animation

% Setup

N_anim=FPS*duration;              % number of frames ( = number of data points)
t_anim=linspace(t0,tf,N_anim);
x_anim=interp1(t_x,x,t_anim);
u_anim=interp1(t_u,u,t_anim)';

pos_anim=x_anim(:,1);                    
v_anim=x_anim(:,2);
I_anim=x_anim(:,3);
T_anim=x_anim(:,4);

% Open figure

AnimFig=figure;
AnimFig.WindowState='maximized';

subplot(4,2,[1 2])   % system
hold on
wire1=plot([-0.08 -0.08],[0.3 -0.1],'Color',hsv2rgb([1 0 T_anim(1)/80]),'LineWidth',1.5);
wire2=plot([-0.04 0],[0.1 -0.1],'Color',hsv2rgb([1 0 T_anim(1)/80]),'LineWidth',1.5);
wire3=plot([0.04 0.08],[0.1 -0.1],'Color',hsv2rgb([1 0 T_anim(1)/80]),'LineWidth',1.5);
pl_mass=rectangle('Position',[pos_anim(1)-0.1 -0.05 0.2 0.1],'FaceColor',[.75 .75 .75]);
wire4=plot([-0.08 -0.04],[-0.1 0.1],'Color',hsv2rgb([1 0 T_anim(1)/80]),'LineWidth',1.5);
wire5=plot([0 0.04],[-0.1 0.1],'Color',hsv2rgb([1 0 T_anim(1)/80]),'LineWidth',1.5);
wire6=plot([0.08 0.08],[-0.1 0.3],'Color',hsv2rgb([1 0 T_anim(1)/80]),'LineWidth',1.5);
rectangle('Position',[pos_anim(end)-0.1 -0.05 0.2 0.1],'FaceColor', [0.75, 0.75, 0.75],'FaceAlpha',0.2,'EdgeColor',[0 0 0 0.2])
grid on
box on
ax_sh=gca;
ax_sh.YGrid='off';
ax_sh.YTickLabel=[];
xlabel('x [m]','Interpreter','latex')
title('Animation')
axis([-abs(max(pos_anim))-0.12 abs(max(pos_anim))+0.12 -0.3 0.3])

subplot(4,2,3)     % position
hold on
pl_1=plot(nan,nan,'b');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Position')
xlim([t0 tf]);
ylim([min(pos_anim) max(pos_anim)])

subplot(4,2,4)     % velocity
hold on
pl_2=plot(nan,nan,'g');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Velocity')
xlim([t0 tf])
ylim([min(v_anim) max(v_anim)])

subplot(4,2,5)     % current
hold on
pl_3=plot(nan,nan,'c');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('I [A]','Interpreter','latex')
title('Current')
xlim([t0 tf])
ylim([min(I_anim) max(I_anim)])

subplot(4,2,6)     % temperature
hold on
pl_4=plot(nan,nan,'m');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('T [$^o$C]','Interpreter','latex')
title('Temperature')
xlim([t0 tf])
ylim([min(T_anim) max(T_anim)])

subplot(4,2,[7 8])     % control u
hold on
pl_u=plot(nan,nan,'r');
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('u [V]','Interpreter','latex')
title('Voltage')
xlim([t0 tf])

pause(2)

% Animation

for ii=2:N_anim

    % compute wire color: up to 80°C, otherwise lampeggia
    if T_anim(ii)>=0 && T_anim(ii)<=80
        value=T_anim(ii)/80; 
    elseif T_anim(ii)<0
        value=0;
    elseif mod(ii,2)==1
        value=1;
    else
        value=0;
    end

    % update plots
    set(pl_mass,'position',[pos_anim(ii)-0.1 -0.05 0.2 0.1])
    set(wire1,'color',hsv2rgb([0 1 value]))
    set(wire2,'color',hsv2rgb([0 1 value]))
    set(wire3,'color',hsv2rgb([0 1 value]))
    set(wire4,'color',hsv2rgb([0 1 value]))
    set(wire5,'color',hsv2rgb([0 1 value]))
    set(wire6,'color',hsv2rgb([0 1 value]))
    set(pl_1,'XData',t_anim(1:ii-1),'YData',pos_anim(1:ii-1))
    set(pl_2,'XData',t_anim(1:ii-1),'YData',v_anim(1:ii-1))
    set(pl_3,'XData',t_anim(1:ii-1),'YData',I_anim(1:ii-1))
    set(pl_4,'XData',t_anim(1:ii-1),'YData',T_anim(1:ii-1))
    set(pl_u,'XData',t_anim(1:ii-1),'YData',u_anim(1:ii-1))

    drawnow
    pause(1/FPS);
    % save GIF
    frame = getframe(AnimFig);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    if ii == 2
        imwrite(A, map, 'animationSimulation.gif', 'gif', 'LoopCount', inf, 'DelayTime', 0.0001);
    else
        imwrite(A, map, 'animationSimulation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.0001);
    end

end
imwrite(A, map, 'animationSimulation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 2);
pause(2)

%% Plots

% plot control in time

figure;
plot(t_u,u,'r')
grid on
xlabel('t [s]','Interpreter','latex')
ylabel('u [V]','Interpreter','latex')
title('Control input (voltage)')
axis tight

% plot state in time

figure;

subplot(2,2,1)
plot(t_x,x(:,1),'b')
grid on
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Position')

subplot(2,2,2)
plot(t_x,x(:,2),'g')
grid on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Velocity')

subplot(2,2,3)
plot(t_x,x(:,3),'c')
grid on
xlabel('t [s]','Interpreter','latex')
ylabel('I [A]','Interpreter','latex')
title('Current')

subplot(2,2,4)
plot(t_x,x(:,4),'m')
grid on
xlabel('t [s]','Interpreter','latex')
ylabel('T [$^o$C]','Interpreter','latex')
title('Temperature')