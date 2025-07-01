function [z_bound]=AssignZbound(x_bound,u_bound,tf_bound,param)

% Takes the values x_bound=[x1_bound, x2_bound, ..., xn_bound] (where n=number of states),
% u_bound=[u1_bound, u2_bound, ..., um_bound] (where m=number of controls) and tf_bound and fills properly
% z_bound with them


% extract parameters

N=param.N;
nx=param.nx;
nu=param.nu;

% creates x_bound and u_bound matrices (each i-th column is the bound for 
% each time instant on the i-th state (or i-th control)

x_bound_mat=zeros(N+1,nx);
for ii=1:nx
    x_bound_mat(:,ii)=x_bound(ii)*ones(N+1,1);
end

u_bound_mat=zeros(N,nu);
for ii=1:nu
    u_bound_mat(:,ii)=u_bound(ii)*ones(N,1);
end

% initialise z_bound (without accounting for tf_bound)

z_bound=zeros(N*(nx+nu)+nx,1); 

% assign u_bound into z_bound

ii=0;
indx_flag=(nx+1)+ii*(nx+nu)+(nu-1);
while indx_flag<length(z_bound)-nx
    z_bound((nx+1)+ii*(nx+nu):(nx+1)+ii*(nx+nu)+(nu-1))=u_bound_mat(ii+1,:);
    indx_flag=(nx+1)+ii*(nx+nu)+(nu-1);
    ii=ii+1;
end

% assign x_bound into z_bound

ii=0;
indx_flag=1+ii*(nx+nu)+(nx-1);
while indx_flag<length(z_bound)
    z_bound(1+ii*(nx+nu):1+ii*(nx+nu)+(nx-1))=x_bound_mat(ii+1,:);
    indx_flag=1+ii*(nx+nu)+(nx-1);
    ii=ii+1;
end

% assign tf_bound into z_bound

z_bound=[z_bound; tf_bound];

end