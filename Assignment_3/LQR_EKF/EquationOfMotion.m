function xp=EquationOfMotion(t,x,u,t_u,param)

% extract parameters

m=param.m;
c=param.c;
k1=param.k1;
k3=param.k3;
alpha=param.alpha;
L0=param.L0;
beta1=param.beta1;
beta2=param.beta2;
RES=param.RES;
Ct=param.Ct;
h=param.h;
Ta=param.Ta;

% ode45 discretises ad cazzum: make sure to take u at the "correct" time
% istant (coherent with x which is being integrated)

u=interp1(t_u,u,t);

% differential equation of motion

xp=zeros(4,1);

xp(1)=x(2);
xp(2)=-k1/m*x(1)-k3/m*x(1)^3-c/m*x(2)+alpha/m*x(3);
xp(3)=-RES/(L0+beta1*x(4)+beta2*x(4)^2)*x(3)+u/(L0+beta1*x(4)+beta2*x(4)^2);
xp(4)=1/Ct*(RES*x(3)^2-h*(x(4)-Ta));

end