% Nicola Visentin, 20/10/2024, Politecnico di Milano

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       ASSIGNMENT 1: SHUTTLE                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Case:   fixed end-state and fixed end-time                            %                                     
% Method: gradient descendent (indirect)                                %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots and animations settings

% show progression of the optimisation (it takes longer)

animation=1;
N_anim=10;    % "smoothness" of the animation (update plots every N_anim iterations)

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
r_control = 0.01/(m_i*g)^2;      % weight for the control action (integral cost)
q_altitude = 0/h_f^2;            % weight for the altitude (integral cost)
q_velocity = 0/1000^2;           % weight for the speed (integral cost)
q_mass = 0/m_i^2;                % weight for the mass (integral cost)
q_angle = 0/(2*pi)^2;            % weight for the inclination (integral cost)

P=[p_altitude 0 0 0; 0 p_velocity 0 0; 0 0 p_mass 0; 0 0 0 p_angle];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
R=r_control;                                                             % L=0.5*(u'*R*u+x'*Q*x)
Q=[q_altitude 0 0 0; 0 q_velocity 0 0; 0 0 q_mass 0; 0 0 0 q_angle];     % L=0.5*(u'*R*u+x'*Q*x)

% ode45 options

options=odeset('RelTol',1e-8,'AbsTol',[1e-8 1e-8 1e-8 1e-8]);

%% Optimisation procedure: gradient descendent method

% set parameters

N=1001;           % number of discretisation intervals

tau_i=1e9;        % speed of control adjustment (exponentially reduces during iterations)
tau_f=1e-2;

ii_max=1001;      % max number of iterations
tol=1e-13;        % tolerance on the optimality condition to be satisfied

% divide the time interval and discretise control u (make a piecewise 
% constant first guess)

t_u=linspace(t0,tf,N);
u=m_i*g*ones(N,1);

% if animation==1, create figure to see the progression of the optimisation

if animation==1

    AnimFig=figure;
    AnimFig.WindowState='maximized';
    
    subplot(3,2,1)
    hold on
    pl_h=plot(nan,nan,'b');
    xlabel('t [s]',Interpreter='latex')
    ylabel('h [m]',Interpreter='latex')
    title('Altitude')
    grid on
    box on
    xlim([t0 tf])
    h_goal=plot(tf,h_f,'bo');
    
    subplot(3,2,2)
    hold on
    pl_v=plot(nan,nan,'g');
    xlabel('t [s]',Interpreter='latex')
    ylabel('v [m/s]',Interpreter='latex')
    title('Velocity')
    grid on
    box on
    xlim([t0 tf])
    v_goal=plot(tf,v_f,'go');

    subplot(3,2,3)
    hold on
    pl_m=plot(nan,nan,'c');
    xlabel('t [s]',Interpreter='latex')
    ylabel('m [kg]',Interpreter='latex')
    title('Mass')
    grid on
    box on
    xlim([t0 tf])
    m_goal=plot(tf,m_f,'co');

    subplot(3,2,4)
    hold on
    pl_gamma=plot(nan,nan,'m');
    xlabel('t [s]',Interpreter='latex')
    ylabel('$\gamma$ [deg]',Interpreter='latex')
    title('Inclination')
    grid on
    box on
    xlim([t0 tf])
    gamma_goal=plot(tf,rad2deg(gamma_f),'mo');

    subplot(3,2,5)
    hold on
    pl_u=plot(nan,nan,'r');
    xlabel('t [s]',Interpreter='latex')
    ylabel('u [N]',Interpreter='latex')
    title('Thrust')
    grid on
    box on
    xlim([t0 tf])

    subplot(3,2,6)
    hold on
    pl_J=semilogy(nan,nan,'.');
    xlabel('iteration',Interpreter='latex')
    ylabel('J',Interpreter='latex')
    title('Cost functional')
    grid on
    box on

end

% iterate until the optimality equation is decently satisfied or the
% iterations limit is reached

ii=-1;
dHdu_norm=tol+eps;

J=zeros(1,ii_max+1);
tic
while ii<ii_max && dHdu_norm>tol

    % set tau parameter

    tau=(tau_i-tau_f)*exp(-5*(ii+1)/ii_max)+tau_f;

    % integrate forward the equation of motion, having u and the 
    % initial condition on x
    
    [t_x,x]=ode45(@(t,x) EquationOfMotion(t,x,u,t_u,param),[t0 tf],x_i,options);
    
    % find final condition for lambda, having u and x at t=tf
    
    x_tf=x(end,:)';
    lambda_f=((x_tf-x_f)'*P)';
    
    % integrate backward the adjoint condition, having x,u and the final
    % condition on lambda
    
    [t_lambda,lambda]=ode45(@(t,lambda) AdjointCondition(t,lambda,u,t_u,x,t_x,param,Q),[tf t0],lambda_f,options);
    lambda=interp1(t_lambda,lambda,t_x);    % this both "reverses" lambda (which originally is defined from tf=t_lambda(1) to t0=t_lambda(end)) and re-samples it consistently with how x is sampled
    
    % check optimality condition
    
    dHdu=Compute_dHdu(R,lambda,u,t_u,x,t_x,param);
    dHdu_norm=norm(dHdu)^2;

    % compute the corresponding cost functional

    J(ii+2) = 0.5*sum(sum(u.^2*R))*tf/length(t_u) + 0.5*(x_tf-x_f).'*P*(x_tf-x_f) + 0.5*sum(sum(x.^2*Q))*tf/length(t_x);
    
    % adjust guess for u (if necessary)
    
    if dHdu_norm>tol
        u_old=u;
        u=AdjustControlGuess(u_old,tau,dHdu,t_x,t_u);
    else
        J(ii+3:end)=[];
    end

    % if animation==1, update the plots every N_anim iterations

    if animation==1 && ((mod(ii,N_anim)==0 || ii==-1) || (mod(ii,N_anim/10)==0 && ii<150))

        set(pl_h,'XData',t_x,'YData',x(:,1))
        set(h_goal,'YData',h_f);
        set(pl_v,'XData',t_x,'YData',x(:,2))
        set(v_goal,'YData',v_f);
        set(pl_m,'XData',t_x,'YData',x(:,3))
        set(m_goal,'YData',m_f);
        set(pl_gamma,'XData',t_x,'YData',rad2deg(x(:,4)))
        set(gamma_goal,'YData',rad2deg(gamma_f));
        set(pl_u,'XData',t_u,'YData',u)
        set(pl_J,'XData',1:ii+2,'YData',J(1:ii+2))
        set(gca,'YScale','log')
        drawnow

        % % save GIF
        % frame = getframe(AnimFig);
        % im = frame2im(frame);
        % [A, map] = rgb2ind(im, 256);
        % if ii == -1
        %     imwrite(A, map, 'animationOptimisation.gif', 'gif', 'LoopCount', inf, 'DelayTime', 0.001);
        % else
        %     imwrite(A, map, 'animationOptimisation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.001);
        % end

    end

    ii=ii+1;

end
ela_time=toc;

if animation==1
    % imwrite(A, map, 'animationOptimisation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 2);
    pause(2)
    close(AnimFig)
end

fprintf('Time required for optimisation: %g s \t (%g/%g iterations)\n',ela_time,ii,ii_max)
fprintf('Initial cost: J=%g\n',J(1))
fprintf('Fianl cost: J=%g\n\n',J(end))

% save optimal control for simulation, if you want

save('optimal_u.mat','u')

%% Plots

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

% plot cost functional through the iterations

figure
semilogy(1:ii+1,J,'.-')
grid on
box on
xlabel('number of iterations',Interpreter='latex')
ylabel('J',Interpreter='latex')
title('Cost functional')
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