function B=dfdu(x,param)

% extract parameters

L0=param.L0;
beta1=param.beta1;
beta2=param.beta2;

% compute dfdu

B=[0 0 1/(L0+beta1*x(4)+beta2*x(4)^2) 0]';

end