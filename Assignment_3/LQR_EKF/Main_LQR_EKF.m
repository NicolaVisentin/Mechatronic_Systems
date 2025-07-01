clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                            ASSIGNMENT 3                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Control:   LQR                                                        %                                                                              
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

%% LQR

% LQR parameters (weights on delta_x and delta_u)

p_delta1 = 0/1^2;
p_delta2 = 0/20^2;
p_delta3 = 0/20^2;
p_delta4 = 0/20^2;
r_delta_control = 10/1^2;
q_delta1 = 100/1^2;
q_delta2 = 1/1^2;
q_delta3 = 1/1^2;
q_delta4 = 0/1^2;

P_lqr=[p_delta1 0 0 0; 0 p_delta2 0 0; 0 0 p_delta3 0; 0 0 0 p_delta4];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
R_lqr=r_delta_control;                                                      % L=0.5*(u'*R*u+x'*Q*x)
Q_lqr=[q_delta1 0 0 0; 0 q_delta2 0 0; 0 0 q_delta3 0; 0 0 0 q_delta4];     % L=0.5*(u'*R*u+x'*Q*x)

% Nominal trajectory

load('opt_data.mat');

u_opt=opt_data.u;
x_opt=opt_data.x; 
t_u=opt_data.t_u;
t_x=opt_data.t_x;

% Re-sample nominal trajectory

N=1001;    % number of time intervals for re-sampling
time=linspace(0,t_x(end),N);

nx=min(size(x_opt));
nu=min(size(u_opt));

xk=interp1(t_x,x_opt,time);
uk=interp1(t_u,u_opt,time)';

% Solve Riccati equation to find matrix PP.
%
%   ! backward integration
%   ! Riccati equation is a matrix equation, but ode45 only works with
%     vectors; we need to "unwrap" all matrices in input to ode45

options=odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,nx^2));
t0=time(1);
tf=time(end);

PP_tf=P_lqr;                % PP(tf) = d^2/dx^2(Phi)|tf = P
PP_tf_vect=P_lqr(1:end)';   % convert initial condition matrix into vector to input ode45

[t_PP,PP_vect]=ode45(@(t,PP) DRE(t,PP,Q_lqr,R_lqr,xk,time,uk,time,param),[tf t0],PP_tf_vect,options);

% Note that PP_vect is actually a (N x nx^2) matrix, where the i-th row
% corresponds to a vector that represents the "unwrapped" Riccati matrix in
% the i-th time instant. Let's flip, re-sample and reshape PP_vect to
% obtain Riccati matrix time history in a (nx x nx x N) 3D matrix where each
% one of the N "slices" corresponds to the nx x nx Riccati matrix at a
% certain time instant

PP_vect=flipud(PP_vect);              % flip PP_vector (was backward integrated)
t_PP=flipud(t_PP);                    % flip t_PP vector (was backward integrated)
PP_vect=interp1(t_PP,PP_vect,time);   % re-sample PP_vect
PP=reshape(PP_vect.',nx,nx,[]);       % reshape PP_vect into PP as said before
PP=permute(PP,[2, 1, 3]);

% Compute the gain matrix K time history ( K is a nx by nu by N 3D matrix; 
% each "slice" nu x nx is the K matrix at one of the N time instants)

B=dfdu(xk,param);
K_lqr=pagemtimes(inv(R_lqr)*B',PP);     % K(t)=R^-1*B(t)*P(t)

%% Extended Kalman filter

% Covariance matrices

Qk=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];   % reliability of the system model
Rk=[1 0; 0 10];                            % reliability of the measurements (sensors)

%% Simulink simulation

% Data for simulink

IC=[0 0 0 20]';                             % nonlinear system initial condition
%IC_obs=IC;                                 % observer initial conditions (exact)
IC_obs=[0.2+IC(1), IC(2), IC(3)+0.001, IC(4)]';  % observer initial conditions (not exact)

w_amp=0.4e-2;    % disturbances amplitude  
npos_amp=5e-2;   % position measurement noise amplitude
ncurr_amp=1e-2;  % current measurement noise amplitude

C=[1 0 0 0; 0 0 1 0];    % output matrix (defines what are the outputs we are measuring)
K_lqr=squeeze(K_lqr)';
time=time';

% Simulate

out=sim('LQR_EKF.slx');

% Extract data

t_sim=squeeze(out.tout);
u=squeeze(out.u);
u_ref=squeeze(out.u_ref);
x_hat=squeeze(out.x_hat)';
x=squeeze(out.x)';
x_ref=squeeze(out.x_ref)';
y=squeeze(out.y);

% Plots

figure    % LQR: actual states vs desired states

subplot(2,2,1)
hold on
plot(t_sim,x(:,1),'b')
plot(t_sim,x_ref(:,1),'r--')
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
plot(t_sim,x_ref(:,2),'r--')
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
plot(t_sim,x_ref(:,3),'r--')
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
plot(t_sim,x_ref(:,4),'r--')
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
plot(t_sim,x(:,1),'b')
plot(t_sim,x_hat(:,1),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [m]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$x$','$\hat{x}$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,2)
hold on
plot(t_sim,x(:,2),'b')
plot(t_sim,x_hat(:,2),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$v$','$\hat{v}$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,3)
hold on
plot(t_sim,x(:,3),'b')
plot(t_sim,x_hat(:,3),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('I [A]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$I$','$\hat{I}$','Interpreter','LaTex','Location','best')
axis tight
hold off

subplot(2,2,4)
hold on
plot(t_sim,x(:,4),'b')
plot(t_sim,x_hat(:,4),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('T [$^o$C]','Interpreter','latex')
title('Observer: observed vs actual')
legend('$T$','$\hat{T}$','Interpreter','LaTex','Location','best')
axis tight
hold off

figure   % observed vs measured vs actual position and current

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

%% Save control for simulation

save('actual_control.mat','u','t_sim');