function dHdu=Compute_dHdu(R,lambda,u,t_u,x,t_x,param)

% extract parameters

alpha=param.alpha;

% ode45 discretises ad cazzum: make sure to take x, u and lambda that are
% discretised at the same way

u=interp1(t_u,u,t_x);

dHdu=R*u+lambda(:,2)./x(:,3)-alpha*lambda(:,3);

end