function A=dfdx(x,u,param)

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

% compute dfdx

A=zeros(4);

A(1,2)=1;
A(2,1)=-(k1+3*k3*x(1)^2)/m;
A(2,2)=-c/m;
A(2,3)=alpha/m;
A(3,3)=-RES/(L0+beta1*x(4)+beta2*x(4)^2);
A(3,4)=((RES*x(3)-u)*(beta1+2*beta2*x(4)))/(L0+beta1*x(4)+beta2*x(4)^2)^2;
A(4,3)=2*RES*x(3)/Ct;
A(4,4)=-h/Ct;

end