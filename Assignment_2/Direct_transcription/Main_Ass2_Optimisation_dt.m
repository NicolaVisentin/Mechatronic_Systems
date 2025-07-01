% Nicola Visentin, 20/10/2024, Politecnico di Milano

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       ASSIGNMENT 2: SHUTTLE                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Case:   fixed end-state and fixed end-time                            %                                                                              
% Method: direct transcription (direct)                                 %
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

x_i=[h_i v_i m_i gamma_i]';
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

% time intervals (discretisation)

N=201;

% solver options

% options=optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'FiniteDifferenceType','central','MaxFunctionEvaluations',1e22,'OptimalityTolerance',1e-10,'StepTolerance',1e-20,'UseParallel',true,'PlotFcn','optimplotfval','Display','iter');
options=optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'FiniteDifferenceType','central','MaxFunctionEvaluations',1e22,'OptimalityTolerance',1e-10,'StepTolerance',1e-20,'UseParallel',true);

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

% define initial guess: z_guess=[x0' u0' x1' u1' ... xN-1' uN-1' xN']'

u_guess=m_i*g*ones(N,nu);             % choose initial guess for the control
z_guess=assign_z0(nx,nu,N,h,x_i,u_guess,xp);

% create a structure to containing parameters and stuff to simplify

param.N=N;
param.nu=nu;
param.nx=nx;
param.xp=xp;
param.x_in=x_i;
param.dfdx=dfdx;
param.dfdu=dfdu;
param.L=L;
param.dLdx=dLdx;
param.dLdu=dLdu;
param.phi=phi;
param.dphidx=dphidx;
param.h=h;

% compute initial cost

J_in=CostAndGrad(z_guess,param);

%% Optimisation

% define objective function 

ObjFun=@(z) CostAndGrad(z,param);

% define bounds on state and control

lb = []; % lower bound
ub = []; % upper bound

% define nonlinear constraints

NLcon=@(z) NonlinConstraintAndGrad(z,param);

% define linear constraints

A = [];
b = [];
Aeq = [];
beq = []; 

% check the supplied gradients

check_options=optimoptions("fmincon",FiniteDifferenceType="central");

[valid_cost,err_gradcost]=checkGradients(@(z) CostAndGrad(z,param), z_guess,check_options,'Display','on');
[valid_con,err_gradcon]=checkGradients(@(z) NonlinConstraintAndGrad(z,param), z_guess,check_options,'IsConstraint',true,'Display','on');

fprintf('CheckGradients on cost function: %g\n',valid_cost)
fprintf('CheckGradients on inequality constraints: %g\n',valid_con(1))
fprintf('CheckGradients on equality constraints: %g\n\n',valid_con(2))

% minimize ObjFun

tic
[z,fval]=fmincon(ObjFun,z_guess,A,b,Aeq,beq,lb,ub,NLcon,options);
ela_time=toc;

fprintf('\nTime required for optimisation: %g s\n',ela_time)
fprintf('Initial cost: J=%g\nFianl cost: J=%g\n\n',J_in,fval)

%% Plot the results

% build time vectors

t_x=0:h:N*h;
t_u=t_x(1:end-1);

% extract states and control from z

x=zeros(N+1,nx); 
u=zeros(N,nu);

for ii=0:N
    x(ii+1,:)=z((1+ii*(nu+nx)):(nx+ii*(nu+nx)));
end    
for ii=0:N-1
    u(ii+1,:)=z((1+nx+ii*(nu+nx)):(nx+nu+ii*(nu+nx)));    
end

% plot optimal control in time

figure
plot(t_u,u,'r')
grid on
box on
xlabel('t [s]',Interpreter='latex')
ylabel('u [N]',Interpreter='latex')
title('Optimal control')
legend('$u^*(t)$','interpreter','latex','location','best')
axis tight

% plot optimal state in time

figure

subplot(2,2,1)
hold on
plot(t_x,x(:,1),'b')
plot(tf,x_f(1),'bo')
grid on
box on
xlabel('t [s]',Interpreter='latex')
ylabel('h [m]',Interpreter='latex')
title('Optimal altitude profile')
legend('$h^*(t)$','target altitude','interpreter','latex','location','best')
axis tight
hold off

subplot(2,2,2)
hold on
plot(t_x,x(:,2),'g')
plot(tf,x_f(2),'go')
grid on
box on
xlabel('t [s]',Interpreter='latex')
ylabel('v [m/s]',Interpreter='latex')
title('Optimal velocity profile')
legend('$v^*(t)$','target speed','interpreter','latex','location','best')
axis tight
hold off

subplot(2,2,3)
hold on
plot(t_x,x(:,3),'c')
plot(tf,x_f(3),'co')
grid on
box on
xlabel('t [s]',Interpreter='latex')
ylabel('m [kg]',Interpreter='latex')
title('Optimal mass profile')
legend('$m^*(t)$','target mass','interpreter','latex','location','best')
axis tight
hold off

subplot(2,2,4)
hold on
plot(t_x,rad2deg(x(:,4)),'m')
plot(tf,rad2deg(x_f(4)),'mo')
grid on
box on
xlabel('t [s]',Interpreter='latex')
ylabel('$\gamma$ [deg]',Interpreter='latex')
title('Optimal angle profile')
legend('$\gamma^*(t)$','target angle','interpreter','latex','location','best')
axis tight
hold off

% save optimal control for simulation, if you want

save('optimal_u.mat','u')