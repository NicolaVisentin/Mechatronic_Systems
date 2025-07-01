function z_guess=assign_z0(nx,nu,N,h,x_in,u_guess,f)

% This function assigns initial guess z_guess for direct methods 
% optimisation
%
%   z_guess=assign_z0(nx,nu,N,x_in,u_guess,f)
%
% Inputs:
%   nx: number of state variables
%   nu: number of control variables
%   N: number of time intervals
%   h: time interval duration
%   x_in: initial condition on the state
%   u_guess: initial guess for the control history
%   f: equation of motion
%
% Output:
%   z_guess: guess for the control and state history, in the form:
%           z_guess=[x0' u0' x1' u1' ... xN-1' uN-1' xN']'
%
% ! ATTENTION: change this function basing on the integration method used (ex.:
% forward euler, crank-nicolson, etc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise z_guess

z_guess=zeros(N*(nx+nu)+nx,1); 

% assign u_guess into z_guess

ii=0;
indx_flag=(nx+1)+ii*(nx+nu)+(nu-1);
while indx_flag<length(z_guess)-nx
    z_guess((nx+1)+ii*(nx+nu):(nx+1)+ii*(nx+nu)+(nu-1))=u_guess(ii+1,:);
    indx_flag=(nx+1)+ii*(nx+nu)+(nu-1);
    ii=ii+1;
end

% compute the consequent x_guess, that depends on u_guess (apart from
% x_in). ATTENTION: forward euler is used here, change it if needed

x_guess=zeros(N+1,nx);
x_guess(1,:)=x_in;

for ii=1:N
    x_ii=x_guess(ii,:)';
    u_ii=u_guess(ii,:)';
    x_guess(ii+1,:)=x_ii+h*f(x_ii,u_ii);
end

% assign x_guess into z_guess

ii=0;
indx_flag=1+ii*(nx+nu)+(nx-1);
while indx_flag<length(z_guess)
    z_guess(1+ii*(nx+nu):1+ii*(nx+nu)+(nx-1))=x_guess(ii+1,:);
    indx_flag=1+ii*(nx+nu)+(nx-1);
    ii=ii+1;
end

end