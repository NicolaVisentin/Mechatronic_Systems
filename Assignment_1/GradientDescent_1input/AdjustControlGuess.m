function u_new=AdjustControlGuess(u_old,tau,dHdu,t_x,t_u)

% resample dHdu such that it is consistent with u

dHdu=interp1(t_x,dHdu,t_u);
dHdu=dHdu';

% update guess

u_new=u_old-tau*dHdu;

end