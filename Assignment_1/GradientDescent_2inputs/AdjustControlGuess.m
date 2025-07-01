function u_new=AdjustControlGuess(u_old,tau,dHdu,t_x,t_u)

% resample dHdu such that it is consistent with u

dHdu=interp1(t_x,dHdu,t_u);

% update guess

u_new(:,1)=u_old(:,1)-tau(1)*dHdu(:,1);
u_new(:,2)=u_old(:,2)-tau(2)*dHdu(:,2);

end