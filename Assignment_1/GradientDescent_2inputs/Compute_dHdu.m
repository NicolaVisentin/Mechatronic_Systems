function dHdu=Compute_dHdu(R,lambda,u,t_u,x,t_x,param)

% extract parameters

alpha=param.alpha;

% ode45 discretises ad cazzum: make sure to take x, u and lambda that are
% discretised at the same way

u=interp1(t_u,u,t_x);

% fill dHdu

dHdu(:,1)=u(:,1)*R(1,1)+u(:,2)*R(2,1)+lambda(:,2)./x(:,3).*cos(u(:,2))-alpha*lambda(:,3)+lambda(:,4).*sin(u(:,2))./x(:,2)./x(:,3);
dHdu(:,2)=u(:,1)*R(1,2)+u(:,2)*R(2,2)-lambda(:,2)./x(:,3).*u(:,1).*sin(u(:,2))+lambda(:,4).*u(:,1).*cos(u(:,2))./x(:,2)./x(:,3);

end