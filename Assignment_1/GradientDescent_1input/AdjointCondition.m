function lambdap=AdjointCondition(t,lambda,u,t_u,x,t_x,param,Q)

% extract parameters

Cd=param.Cd;
Cl=param.Cl;
g=param.g;
rho=param.rho;
S=param.S;

% ode45 discretises ad cazzum: make sure to take x and u at the "correct" 
% time istant (coherent with lambda which is being integrated)

x=interp1(t_x,x,t);
u=interp1(t_u,u,t);

% adjoint differential equation

df1dx1 = 0;
df1dx2 = sin(x(4));
df1dx3 = 0;
df1dx4 = x(2)*cos(x(4));
df2dx1 = 0;
df2dx2 = -rho*S*Cd*x(2)/x(3);
df2dx3 = rho*S/2*Cd*x(2)^2/(x(3))^2-u/x(3)^2;
df2dx4 = -g*cos(x(4));
df3dx1 = 0;
df3dx2 = 0;
df3dx3 = 0;
df3dx4 = 0;
df4dx1 = 0;
df4dx2 = rho*S/2*Cl/x(3)+g*cos(x(4))/x(2)^2;
df4dx3 = -rho*S/2*Cl*x(2)/x(3)^2;
df4dx4 = g*sin(x(4))/x(2);

dfdx=[df1dx1 df1dx2 df1dx3 df1dx4; 
      df2dx1 df2dx2 df2dx3 df2dx4;
      df3dx1 df3dx2 df3dx3 df3dx4
      df4dx1 df4dx2 df4dx3 df4dx4];

dLdx=[x(1)*Q(1,1) x(2)*Q(2,2) x(3)*Q(3,3) x(4)*Q(4,4)];

lambdap=-(dLdx+lambda'*dfdx)';

end