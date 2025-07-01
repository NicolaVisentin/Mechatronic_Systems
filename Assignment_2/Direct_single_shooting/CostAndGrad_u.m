function [J,DJ] = CostAndGrad_u(u_vect,param)

% Gives the cost function J and, if requested, also its gradient

% extract parameters

N = param.N;
nu = param.nu;
nx = param.nx;
h = param.h;
L = param.L;
phi = param.phi;
xp=param.xp;
x_in=param.x_in;

% compute, for each u(t), the corresponding x(t) with forward Euler

x_vect=zeros(nx,N+1);
x_vect(:,1)=x_in;
for ii=1:N
    x_vect(:,ii+1)=x_vect(:,ii)+h*xp(x_vect(:,ii),u_vect(:,ii));
end
    
% compute the value of cost function J

J = 0;
for ii = 1:N
    J = J + h*L(x_vect(:,ii),u_vect(:,ii));
end
J = J + phi(x_vect(:,end));

% compute the value of the gradient of J

if nargout > 1

    dLdx = param.dLdx;
    dLdu = param.dLdu;
    dphidx = param.dphidx;

    % compute the (nx,nu,N+1,N) 4D matrix that collects all sensitivities 
    % dxk/dui (see "sensitivity_dxdu" function for details)
    dxdu=sensitivity_dxdu(x_vect,u_vect,param);
    
    % compute the gradient of J
    DJ=zeros(nu,N);
    for ii = 1:N
        
        DJ(:,ii)=h*dLdu(x_vect(:,ii),u_vect(:,ii));
        for k=ii+1:N
            DJ(:,ii)=DJ(:,ii)+h*dLdx(x_vect(:,k),u_vect(:,k))*dxdu(:,:,k,ii);
        end
        DJ(:,ii)=DJ(:,ii)+dphidx(x_vect(:,N+1),u_vect(:,k))*dxdu(:,:,N+1,ii);

    end

    % vectorise the gradient
    DJ=DJ(:);
    
end

end