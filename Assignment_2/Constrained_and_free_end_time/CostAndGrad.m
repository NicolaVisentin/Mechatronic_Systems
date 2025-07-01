function [J,DJ] = CostAndGrad(z,param)

% Gives the cost function J and, if requested, also its gradient

% extract parameters

N = param.N;
nx = param.nx;
nu = param.nu;
L = param.L;
phi = param.phi;

% extract final time

tf=z(end);
z(end)=[];   % remove tf

% extract states and control from z

x = zeros(nx,N+1); u = zeros(nu,N);
for ii = 0:N
    x(:,ii+1) = z((1 + ii*(nu + nx)):(nx + ii*(nu + nx)));
end    
for ii = 0:N-1
    u(:,ii+1) = z((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)));    
end
    
% compute the value of cost function J

h=tf/N;

J = 0;
for ii = 1:N
    J = J + h*L(x(:,ii),u(:,ii));
end
J = J + phi(x(:,end),tf);

% compute the value of the gradient of J

if nargout > 1

    dLdx = param.dLdx;
    dLdu = param.dLdu;
    dphidx = param.dphidx;
    dphidtf = param.dphidtf;
    L = param.L;
    
    DJ = zeros(size(z));   % without counting tf contribution
    DJ_tf=0;
    for ii = 0:N-1
        DJ((1 + ii*(nu + nx)):(nx + ii*(nu + nx)),1) = h*dLdx(x(:,ii + 1),u(:,ii + 1));
        DJ((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)),1) = h*dLdu(x(:,ii + 1),u(:,ii + 1));   
        DJ_tf = DJ_tf + L(x(:,ii+1),u(:,ii+1));  % update contribution for tf (DJ(end))
    end
    DJ(end - nx + 1:end,1) = dphidx(x(:,end));
    
    DJ=[DJ; dphidtf+1/N*DJ_tf];      % now add tf contribution
    
end

end