function dxdu=sensitivity_dxdu(x_vect,u_vect,param)

% compute sensitivity dx/du of state vector wrt input vector. It's a 4D
% matrix (a,b,c,d) organised like this:
%   - each (:,:,c,d) 2D "face" represents the derivative of state vector at
%     time instant tc wrt input vector at time instant td
%   - each (:,:,c,d) 2D "face" is an (a,b) matrix where a-th rows are
%     associated to states and b-th rows to controls
%   - of course, since each xc state is influenced by all PREVIOUS ud
%     controls (c-1, c-2, c-3, ...), we expect all (:,:,c,d) "2D faces" to
%     null matrices if c<=d

% extract parameters

N=param.N;
nx=param.nx;
nu=param.nu;
dfdx=param.dfdx;
dfdu=param.dfdu;
h=param.h;

% compute the 4D tensor that collects all sensitivities dxc/dud

dxdu=zeros(nx,nu,N+1,N);
for d=1:N
    
    % fill the "easy" "faces" where x=x_i+1 and u=u_i (derivative of the
    % (i+1)-th state vector wrt the previous i-th control vector is
    % simply h*dfdu)
    dxdu(:,:,d+1,d)=h*dfdu(x_vect(:,d+1),u_vect(:,d));

    % fill the other non-null faces (those where c>d)
    for c= d+2 : N+1
        dxdu(:,:,c,d) = (eye(nx)+h*dfdx(x_vect(:,c-1),u_vect(:,c-1))) * dxdu(:,:,c-1,d);
    end

end