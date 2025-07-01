% Nicola Visentin, 03/11/2024, Politecnico di Milano

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       ASSIGNMENT 2: SHUTTLE                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Case:   fixed end-state and fixed end-time; unconstrained             %                                                     %
% Method: direct single shooting (direct)                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data and parameters definition

Cd=0.5;            % drag coefficient [-]
Cl=0.2;            % lift coefficient [-]
alpha=0.00005;     % fuel consumption constant [(kg/s)/N]
g=9.81;            % gravity [m/s^2]
S=2;               % areodynamic surface (hyp.: S=2)
rho=1;             % air density (hyp.: rho=1=const)

param.Cd=Cd;
param.Cl=Cl;
param.alpha=alpha;
param.g=g;
param.S=S;
param.rho=rho;

% optimisation problem data

t0=0;           % initial time [s]
tf=110;         % final time   [s]

h_i=100;        % initial altitude [m]
h_f=50000;      % final altitude [m]

v_i=100;        % initial velocity [m/s]
v_f=1000;       % final velocity [m/s]

m_i=20000;      % initial mass (shuttle + fuel) [kg]
m_f=10000;      % final mass (shuttle + fuel) [kg]

gamma_i=pi/3;   % initial angle [rad]
gamma_f=pi/6;   % final angle [rad]

x_in=[h_i v_i m_i gamma_i]';
x_f=[h_f v_f m_f gamma_f]';

% weights for the control and the state

p_altitude = 5e2/h_f^2;          % weight for the final altitude
p_velocity = 1e3/1000^2;         % weight for the final velocity
p_mass = 1e-1/m_i^2;             % weight for the final mass
p_angle = 1/(2*pi)^2;            % weight for the final inclination angle
r_control = 1e-2/(m_i*g)^2;      % weight for the control action (integral cost)
q_altitude = 0.1/h_f^2;            % weight for the altitude (integral cost)
q_velocity = 0.1/1000^2;           % weight for the speed (integral cost)
q_mass = 0/m_i^2;                % weight for the mass (integral cost)
q_angle = 0/(2*pi)^2;            % weight for the inclination (integral cost)

P=[p_altitude 0 0 0; 0 p_velocity 0 0; 0 0 p_mass 0; 0 0 0 p_angle];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
R=r_control;                                                             % L=0.5*(u'*R*u+x'*Q*x)
Q=[q_altitude 0 0 0; 0 q_velocity 0 0; 0 0 q_mass 0; 0 0 0 q_angle];     % L=0.5*(u'*R*u+x'*Q*x)

% time intervals

N=201;

% solver options

options=optimoptions('fminunc','FiniteDifferenceType','central','MaxFunctionEvaluations',1e20,'OptimalityTolerance',1e-20,'StepTolerance',1e-20,'PlotFcn','optimplotfval','Display','iter');
%options=optimoptions('fminunc','SpecifyObjectiveGradient',true,'PlotFcn','optimplotfval','Display','iter');

%% Define required equations

% equation of motion and its Jacobians

xp=@(x,u) [ x(2)*sin(x(4));
            -rho*S/2*Cd*x(2)^2/x(3)-g*sin(x(4))+u/(x(3));
            -alpha*u;
            rho*S/2*Cl*x(2)/x(3)-g*cos(x(4))/x(2);
            ];

dfdx = @(x,u) dfdx(x,u,param);
dfdu = @(x,u) [ 0 ; 1/x(3); -alpha; 0 ];

% cost functions and their Jacobians

L = @(x,u) 0.5*x'*Q*x+0.5*u'*R*u;
phi = @(x) 0.5*(x-x_f)'*P*(x-x_f);

dLdx = @(x,u) x'*Q;
dLdu = @(x,u) u'*R;
dphidx =@(x) (x-x_f)'*P;

%% Setup for the method

nx=4;           % Number of states
nu=1;           % number of inputs
h=tf/N;         % temporal discretization

% define initial guess [u0' u1' ... uN-2' uN-1']

u_guess=m_i*g*ones(1,N*nu);

% create a structure to containing parameters and stuff to simplify

param.N=N;
param.nu=nu;
param.nx=nx;
param.xp=xp;
param.x_in=x_in;
param.dfdx=dfdx;
param.dfdu=dfdu;
param.L=L;
param.dLdx=dLdx;
param.dLdu=dLdu;
param.phi=phi;
param.dphidx=dphidx;
param.h=h;

%% Optimisation

% check the supplied gradients

u_check=reshape(u_guess,nu,N)';
z_check=assign_z0(nx,nu,N,h,x_in,u_check,xp);   % ! built with forward Euler
check_options=optimoptions("fminunc",FiniteDifferenceType="central");
[valid_cost,err_gradcost]=checkGradients(@(z) CostAndGrad_z(z,param), z_check,check_options,'Display','on');

fprintf('CheckGradients on cost function: %g\n',valid_cost)

% minimize cost J

J_in=CostAndGrad_z(z_check,param);

tic
[u_opt,fval] = fminunc(@(u) CostAndGrad_u(u,param), u_guess, options);
ela_time=toc;

fprintf('\nTime required for optimisation: %g s\n',ela_time)
fprintf('Initial cost: J=%g\nFianl cost: J=%g\n\n',J_in,fval)

%% Plot the results

% build time vector

t=0:h:N*h;

% compute corresponding x out of optimal u, using forward Euler

u_opt=reshape(u_opt,nu,N);
x_opt=zeros(nx,N+1);
x_opt(:,1)=x_in;
for ii=1:N
    x_opt(:,ii+1)=x_opt(:,ii)+h*xp(x_opt(:,ii),u_opt(:,ii));
end

u_opt=u_opt';
x_opt=x_opt';

% plot optimal state and control in time

figure

subplot(2,2,1)
hold on
plot(t,x_opt(:,1),'b')
plot(t(end),x_f(1),'bo')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('h [m]','Interpreter','latex')
title('Optimal altitude profile')
legend('h*(t)','target','interpreter','latex','location','best')
axis tight
hold off

subplot(2,2,2)
hold on
plot(t,x_opt(:,2),'g')
plot(t(end),x_f(2),'go')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('v [m/s]','Interpreter','latex')
title('Optimal velocity profile')
legend('v*(t)','target','interpreter','latex','location','best')
axis tight
hold off

subplot(2,2,3)
hold on
plot(t,x_opt(:,3),'c')
plot(t(end),x_f(3),'co')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('m [kg]','Interpreter','latex')
title('Optimal mass profile')
legend('m*(t)','target','interpreter','latex','location','best')
axis tight
hold off

subplot(2,2,4)
hold on
plot(t,rad2deg(x_opt(:,4)),'m')
plot(t(end),rad2deg(x_f(4)),'mo')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\gamma$ [deg]','Interpreter','latex')
title('Optimal angle profile')
legend('$\gamma$*(t)','target','interpreter','latex','location','best')
axis tight
hold off

figure

plot(t(1:end-1),u_opt,'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('u [N]','Interpreter','latex')
title('Optimal control')
legend('u*(t)','interpreter','latex','Location','best')
axis tight

% save optimal control for simulation, if you want

u=u_opt;
save('optimal_u.mat','u')