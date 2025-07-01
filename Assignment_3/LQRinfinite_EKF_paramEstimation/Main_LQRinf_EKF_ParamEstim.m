clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                            ASSIGNMENT 3                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Control:   Steady-state LQR                                           %                                                                              
% Observer:  Extended Kalman filter                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data and parameters

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

%% Infinite time LQR control definition

% Infinite time LQR parameters (weights on delta_x and delta_u; final target)

r_delta_control = 10/1^2;  % control
q_delta1 = 100/1^2;        % position  
q_delta2 = 1/1^2;          % speed
q_delta3 = 1/1^2;          % current
q_delta4 = 0/1^2;          % temperature

R_lqr=r_delta_control;                                                    % L=0.5*(u'*R*u+x'*Q*x)
Q_lqr=[q_delta1 0 0 0; 0 q_delta2 0 0; 0 0 q_delta3 0; 0 0 0 q_delta4];   % L=0.5*(u'*R*u+x'*Q*x)

x1_f=0.5;
x2_f=0;
x3_f=(k1*x1_f+k3*x1_f^3)/alpha;
x4_f=(RES*x3_f^2+h*Ta)/h;

x_f=[x1_f x2_f x3_f x4_f]';
u_f=RES*x3_f;

nx=size(Q_lqr,1);
nu=size(R_lqr,1);

% A, B, C and D matrices (uncontrolled system)

A=dfdx(x_f,u_f,param);
B=dfdu(x_f,param);
C=[1 0 0 0; 0 0 1 0];
D=0;

poles_un=eig(A);

% Design of the LQR using matlab function

[K_lqr,PP,poles_con]=lqr(A,B,Q_lqr,R_lqr);

% Plot poles (around eq position)

figure
hold on
plot(real(poles_un),imag(poles_un),'bx','LineWidth',1.5,'MarkerSize',10)
plot(real(poles_con),imag(poles_con),'rx','LineWidth',1.5,'MarkerSize',5)
grid on
box on
xlabel('Re','Interpreter','latex')
ylabel('Im','Interpreter','latex')
title('Poles near final position')
legend('poles unc','poles con')
hold off

%% Extended Kalman filter for parameter estimation

% Covariance matrices

Qk=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];   % reliability of the system model
Rk=[10 0; 0 100];                          % reliability of the measurements (sensors)

% Artificial noise covariance on the parameter

Qp=1e8;

% Augmented system

Ca=[1 0 0 0 0; 0 0 1 0 0];   % output matrix
Ra=Rk;                       % output covariance
Qa=[Qk zeros(size(Qk,1),1); zeros(1,size(Qk,2)) Qp];

%% Simulink simulation

% Data for simulink

tf=5;        % simulation time

IC=[0 0 0 20]';                                 % nonlinear system initial condition
%IC_obs=[IC; 3];                                % observer initial conditions (exact, except k1, obv)
IC_obs=[0.02+IC(1), IC(2), IC(3)+0.01, IC(4), 3]';   % observer initial conditions (not exact)

w_amp=0.4e-2;    % disturbances amplitude  
npos_amp=3e-3;   % position measurement noise amplitude
ncurr_amp=2e-2;  % current measurement noise amplitude

% Simulate

out=sim('LQRinf_EKF_ParamEstim.slx');

% Extract data

t_sim=squeeze(out.tout);
u=squeeze(out.u);
u_opt=squeeze(out.u_ref);
x_hat=squeeze(out.x_hat)';
x=squeeze(out.x)';
x_opt=squeeze(out.x_ref)';
y=squeeze(out.y);

% Plots

figure    % LQR: actual states vs desired states

subplot(2,2,1)
hold on
plot(t_sim,x(:,1),'b')
plot(t_sim,ones(1,length(t_sim))*x_opt(1),'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Actual vs reference')
legend('$x$','$x^*$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,2)
hold on
plot(t_sim,x(:,2),'b')
plot(t_sim,ones(1,length(t_sim))*x_opt(2),'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Actual vs reference')
legend('$v$','$v^*$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,3)
hold on
plot(t_sim,x(:,3),'b')
plot(t_sim,ones(1,length(t_sim))*x_opt(3),'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('I [A]','Interpreter','latex')
title('Actual vs reference')
legend('$I$','$I^*$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,4)
hold on
plot(t_sim,x(:,4),'b')
plot(t_sim,ones(1,length(t_sim))*x_opt(4),'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('T [$^o$C]','Interpreter','latex')
title('Actual vs reference')
legend('$T$','$T^*$','Interpreter','LaTex','Location','best')
axis tight
hold off

figure    % OBSERVER: observed states vs real states

subplot(2,2,1)
hold on
plot(t_sim,x_hat(:,1),'r')
plot(t_sim,x(:,1),'b')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$\hat{x}$','$x$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,2)
hold on
plot(t_sim,x_hat(:,2),'r')
plot(t_sim,x(:,2),'b')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$\hat{v}$','$v$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,3)
hold on
plot(t_sim,x_hat(:,3),'r')
plot(t_sim,x(:,3),'b')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('I [A]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$\hat{I}$','$I$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,4)
hold on
plot(t_sim,x_hat(:,4),'r')
plot(t_sim,x(:,4),'b')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('T [$^o$C]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$\hat{T}$','$T$','Interpreter','LaTex','Location','best')
axis tight
hold off

figure   % observed vs measured vs actual position and  current

subplot(211)
hold on
plot(t_sim,y(:,1),'color',[0 1 0 0.3])
plot(t_sim,x(:,1),'b')
plot(t_sim,x_hat(:,1),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Observer: estimation vs measured vs actual position')
legend('$x_{meas}$','$x$','$\hat{x}$','Interpreter','LaTex','Location','best','fontsize',12)
axis tight
hold off

subplot(212)
hold on
plot(t_sim,y(:,2),'color',[0 1 0 0.3])
plot(t_sim,x(:,3),'b')
plot(t_sim,x_hat(:,3),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('I [A]','Interpreter','latex')
title('Observer: estimation vs measured vs actual current')
legend('$I_{meas}$','$I$','$\hat{I}$','Interpreter','LaTex','Location','best','fontsize',12)
axis tight
hold off

figure    % LQR output: u and u*
hold on
plot(t_sim,u,'b')
plot(t_sim,u_opt*ones(length(t_sim),1),'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$u$ [V]','Interpreter','latex')
legend('$u$','$u^*$','Interpreter','latex','Location','best')
axis tight
title('LQR: control voltage')
hold off

% Parameter estimation

figure
hold on
plot(t_sim,x_hat(:,5),'b')
plot(t_sim,k1*ones(length(t_sim),1),'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$k_1$ [N/m]','Interpreter','latex')
axis tight
title('Parameter estimation')
legend('$\hat{k}_1$','$k_1$','Interpreter','latex','Location','best')
hold off

%% Save control for simulation

save('actual_control.mat','u','t_sim');