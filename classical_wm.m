function [u_sol, y_sol, tau_w, T_sol, qw] = classical_wm(y1,u1,rho1,mu1,T1,Tw,case_data,Pr_lam,Pr_t,visc_law,visc_arg)

% get variables from case_data
gamma = case_data.gamma;
R = case_data.R;
Cp = gamma*R/(gamma -1);

tau_w = mu1*u1/y1; % initial guess
qw = mu1*(Cp*T1 + 0.5*u1*u1 - Cp*Tw)/y1; % initial guess
kappa = 0.41;
A_plus = 17;
relax = 0.3;
tol = 1.0e-8;
step_ctl_fax = 0.1;
max_iter = 500;
u_end_prev = nan; tau_w_prev = nan; qw_prev = nan; % these will be set at end of 1st iter

Ny = 3000;
y_sol = zeros(Ny,1);
u_sol = zeros(Ny,1);
T_sol = zeros(Ny,1);

% outer iterations for updating tau_w and qw
for iter = 1:max_iter
    assert(~isnan(tau_w))
    assert(imag(tau_w)==0)
    assert(tau_w>0)
    assert(~isnan(qw))
    assert(imag(qw)==0)

    % assume constant pressure in the inner layer to compute a wall density
    rhow = rho1*T1/Tw;
    if strcmp(visc_law,'power law')
        muw = mu1*power(Tw/T1, visc_arg);
    elseif strcmp(visc_law,'sutherland')
        S = visc_arg;
        C1 = mu1*(T1+S)/T1^(3/2);
        muw = C1*Tw.^(3/2)./(Tw+S);
    else
        assert(0) % unrecognized visc_law input argument
    end
    utau = sqrt(tau_w/rhow);
    l_visc = muw/(rhow*utau); % viscous unit
    y_sol(1)                   = 0.0; % integration coordinate
    int1        = 0.0;
    int2        = 0.0;
    int3        = 0.0;
    int4        = 0.0;
    u_sol(1)           = 0.0; % integrated value
    T_sol(1)           = Tw;
    dy           = 0.1*l_visc; % step size in viscous units
    iter_sec = 1; % integration index
    while ( y_sol(iter_sec) < y1 )
        dy                     = min(dy, y1-y_sol(iter_sec)); % constant step size except last step
        y_half    = y_sol(iter_sec) + 0.5*dy; % current cell center
        mu_t      = kappa*rhow*utau*sqrt(Tw/T_sol(iter_sec))*y_half*(1.0- exp(-y_half/l_visc/A_plus)).^2;
        if strcmp(visc_law,'power law')
            mu_lam = mu1*power(T_sol(iter_sec)/T1, visc_arg);
        elseif strcmp(visc_law,'sutherland')
            S = visc_arg;
            C1 = mu1*(T1+S)/T1^(3/2);
            mu_lam = C1*T_sol(iter_sec).^(3/2)./(T_sol(iter_sec)+S);
        else
            assert(0) % unrecognized visc_law input argument
        end
        mu_tot    = mu_lam + mu_t;
        k_lam     = Cp*mu_lam/Pr_lam;
        k_t       = Cp*mu_t/Pr_t;
        k_tot     = k_lam + k_t;
        unew      = u_sol(iter_sec) + dy*tau_w/mu_tot;
        Tnew      = T_sol(iter_sec) + dy/k_tot*(qw - 0.5*tau_w*u_sol(iter_sec) - 0.5*tau_w*unew);
        % update the integrals ..
        mu_prime  = Pr_t*mu_lam/Pr_lam;
        int1 = int1 + dy/(mu_prime + mu_t);
        int2 = int2 + 0.5*(u_sol(iter_sec) + unew)*dy/(mu_prime + mu_t);
        int3 = int3 + dy/(mu_lam + mu_t);
        int4 = int4 + 0.5*mu_t*dy/(mu_tot*mu_tot)/sqrt(tau_w);
        
        % update the solutions and y
        dudy = (unew-u_sol(iter_sec))/dy;
        dTdy = (Tnew-T_sol(iter_sec))/dy;
        y_sol(iter_sec+1) = y_sol(iter_sec) + dy;
        u_sol(iter_sec+1) = unew;
        T_sol(iter_sec+1) = Tnew;
        iter_sec = iter_sec+1;
        if iter_sec > Ny
            assert(0);
        end

        % adaptive y step
        dy = max(dy,step_ctl_fax*u_sol(iter_sec)/abs(dudy));

    end
    
    y_sol = y_sol(1:iter_sec);
    u_sol = u_sol(1:iter_sec);
    T_sol = T_sol(1:iter_sec);

    assert( abs(y_sol(iter_sec)-y1)/y1 < 1.0e-12);
    % check for convergence ..
    if ( ( (abs(u1) < 1.0e-12) || ( abs(u_sol(iter_sec)-u1)/u1 < tol)))
        break;
    end
    
    % update the estimates for tau_w
    J                   = int3 - sqrt(tau_w)*int4;
    dtau                = (u1-u_sol(iter_sec))/J;
    tau_w           = tau_w + relax*dtau;
    % For isothermal flows:
    dhstar = Cp*(T1-Tw);
    qw = relax* ((dhstar/Pr_t + tau_w*int2)/int1) + (1.0-relax)*qw;
end

if ( iter == max_iter)
    fprintf('reached max_iter!!!!\n')
end
qw = -qw;    
end