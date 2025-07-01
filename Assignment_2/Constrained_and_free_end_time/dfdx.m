function A=dfdx(x,u,param)

% extract parameters

Cd=param.Cd;
Cl=param.Cl;
g=param.g;
rho=param.rho;
S=param.S;

% partial derivatives

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

% assemble

A=[df1dx1 df1dx2 df1dx3 df1dx4; 
      df2dx1 df2dx2 df2dx3 df2dx4;
      df3dx1 df3dx2 df3dx3 df3dx4
      df4dx1 df4dx2 df4dx3 df4dx4];

end