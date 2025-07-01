function xp=EquationOfMotion(t,x,u,t_u,param)

% extract parameters

Cd=param.Cd;
Cl=param.Cl;
alpha=param.alpha;
g=param.g;
rho=param.rho;
S=param.S;

% ode45 discretises ad cazzum: make sure to take u at the "correct" time
% istant (coherent with x which is being integrated)

u=interp1(t_u,u,t);

% differential equation of motion

xp=zeros(4,1);

xp(1)=x(2)*sin(x(4));
xp(2)=u(1)/x(3)*cos(u(2))-rho*S*Cd*x(2)^2/(2*x(3))-g*sin(x(4));
xp(3)=-alpha*u(1);
xp(4)=rho*S*Cl*x(2)/(2*x(3))-g*cos(x(4))/x(2)+u(1)*sin(u(2))/x(2)/x(3);

end