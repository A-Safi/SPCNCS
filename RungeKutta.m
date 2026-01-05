function nstate = RungeKutta(fnc,state,h,u)

[fx,gx] = fnc(state);
k1 = fx+gx*u;

[fx,gx] = fnc(state+h/2*k1);
k2 = fx+gx*u;

[fx,gx] = fnc(state+h/2*k2);
k3 = fx+gx*u;

[fx,gx] = fnc(state+h*k3);
k4 = fx+gx*u;

nstate=state+h/6*(k1+2*k2+2*k3+k4);

end

