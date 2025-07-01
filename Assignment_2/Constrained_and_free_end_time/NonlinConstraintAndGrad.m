function [C,Ceq,Dc,DCeq] = NonlinConstraintAndGrad(z,param)

% computes the value of nonlinear inequality constraint C(z)<=0 and nonlinear
% equality constraint Ceq(z)=0 for a given z, and, if required, also the 
% value of their gradients
%
%   [C,Ceq,Dc,DCeq] = NonlinConstraintAndGrad(z,param)
%
% Outputs
%   C:       value of C in nonlinear inequality constraint C(z)<0
%   Ceq:     value of Ceq in equality nonlinear constraint Ceq(z)=0
%   Dc:      value of the gradient of the inequality nonlinear constraint Dc(z)
%   DCeq:    value of the gradient of the equality nonlinear constraint Dcon(z)
 
% extract parameters

N = param.N;
nx = param.nx;
nu = param.nu;
xp = param.xp;
x_in = param.x_in;

tf=z(end);
h=tf/N;

% extract states and control from z

x = zeros(nx,N+1); 
u = zeros(nu,N);
for ii = 0:N
    x(:,ii+1) = z((1 + ii*(nu + nx)):(nx + ii*(nu + nx)));
end    
for ii = 0:N-1
    u(:,ii+1) = z((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)));    
end

% constraint functions

C = []; % here is room for inequality constraint
Ceq = zeros((N+1)*nx,1);
Ceq(1:nx) = x_in - x(:,1) ; % initial condition constraint

for ii = 1:N
    Ceq((1:nx) + ii*nx) = x(:,ii) + h*xp(x(:,ii),u(:,ii)) - x(:,ii+1);    
end

% gradients of the constraints: it's a ( (N+1)*nx, N*(nu*nx)+nx+1 ) matrix

dfdx = param.dfdx;
dfdu = param.dfdu;

if nargout > 2

    Dc = []; % here is room for inequality constraint
    DCeq = zeros((N+1)*nx, N*(nu+nx) + nx + 1);
    DCeq(1:nx,1:nx) = - eye(nx); 
    for ii = 1:N
        DCeq((1 + nx +(ii - 1)*nx):((nx + ii*nx)),(1 +(ii - 1)*(nx+nu)):((ii - 1)*(nx+nu) + nx)) = eye(nx) + h*dfdx(x(:,ii),u(:,ii));
        DCeq((1 + nx +(ii - 1)*nx):((nx + ii*nx)),(1 +(ii)*(nx+nu)):((ii)*(nx+nu) + nx)) = - eye(nx);
        DCeq((1 + nx +(ii - 1)*nx):((nx + ii*nx)),(1 +(ii - 1)*(nx+nu) + nx):((ii - 1)*(nx+nu) + nx + nu)) = h*dfdu(x(:,ii),u(:,ii));
        DCeq((1 + nx +(ii - 1)*nx):((nx + ii*nx)),end) = 1/N*xp(x(:,ii),u(:,ii));
    end
    DCeq = DCeq.';   % transpose it because ode works with time along the columns, not rows
    
end

end