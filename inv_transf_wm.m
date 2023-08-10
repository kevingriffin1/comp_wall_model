function [u_sol, y_sol, tau_w, T_sol, qw] = inv_transf_wm(y1,u1,rho1,mu1,T1,Tw,case_data,r,visc_law,visc_arg,s_Pr_input)

% solver parameters
max_iter = 300;
max_iter_sec = 300;
relax_u = 1;
relax_T = 1;
tol_u = 1.0e-4;
tol_T = 1.0e-7;
step_ctl_fax = 0.1;

gamma = case_data.gamma;
R = case_data.R;
Cp = gamma*R/(gamma -1); % perfect gas
if isfield(case_data,'Pr')
    Pr = case_data.Pr;
else
    Pr = 0.70;
end
Me = case_data.Me;
Te = case_data.Te;
Ue = case_data.Ue;

assert(isfloat(s_Pr_input))
Tr = Te*(1+r*(gamma-1)/2*Me^2);
rhow = rho1*T1./Tw;
if strcmp('power law',visc_law)
    muw = mu1*(Tw/T1)^visc_arg;
    mue = mu1*(Te/T1)^visc_arg;
elseif strcmp('sutherland',visc_law)
    % Sutherland's viscosity law
    %mu/muref = (T/Tref).^(3/2).*(Tref+S)./(T+S)  %T0=273.15; mu0=1.716e-5;
    S = visc_arg;
    C1 = mu1*(T1+S)/T1^(3/2);
    muw = C1*Tw.^(3/2)./(Tw+S);
    mue = C1*Te.^(3/2)./(Te+S);
else
    assert(0) % unrecognized visc_law input argument
end
rhoe = rho1*T1./Te;
tau_w = muw*u1/y1; % initial guess
utau = sqrt(tau_w/rhow); % initial guess
U_end_prev = nan;
utau_prev = utau;

Ny = 3000;
y_sol = zeros(Ny,1);
u_sol = zeros(Ny,1);
T_sol = zeros(Ny,1);

for iter=1:max_iter % tau_w iteration
    % initialize the data for the point behind in space (used for
    % derivatives)
    T_sol(1) = Tw;
    u_sol(1) = 0;
    rho_m1 = rhow;
    mu_m1 = muw;
    y_sl_m1 = 0;
    i_y = 1;

    l_visc = muw/rhow/utau;
    dy  = l_visc*0.2;
    
    while y_sol(i_y) < y1
        dy = min(dy,y1-y_sol(i_y));
        y_half = y_sol(i_y) + 0.5*dy;
        i_y = i_y + 1;
        y_sol(i_y) = y_sol(i_y-1) + dy;
        T_sol(i_y) = T_sol(i_y-1); % initial guess
        T_prev = T_sol(i_y);
        for iter_sec=1:max_iter_sec % temperature iteration
            rho = Tw/T_sol(i_y)*rhow;
            if strcmp('power law',visc_law)
                mu = (T_sol(i_y)/Tw)^(visc_arg)*muw;
            elseif strcmp('sutherland',visc_law)
                mu = C1*T_sol(i_y).^(3/2)./(T_sol(i_y)+S);
            else
                assert(0) % unrecognized visc_law input argument
            end            
            rho_half = 0.5*(rho+rho_m1);
            mu_half = 0.5*(mu+mu_m1);
            y_sl = y_sol(i_y)*utau.*sqrt(rho*rhow)./mu;
            dysl = y_sl - y_sl_m1;
            drho_dy = (rho-rho_m1)/dy;
            dmu_dy = (mu-mu_m1)/dy;
            f = muw/mu_half - sqrt(rho_half/rhow).*(1+0.5*drho_dy.*y_half./rho_half - dmu_dy.*y_half./mu_half);
            kappa = 0.41; 
            A_plus = 17;
            y_sl_half = (y_sl + y_sl_m1)/2;
            mu_t      = (kappa*y_sl_half*(1.0- exp(-y_sl_half/A_plus)).^2);
            dU_incomp_dysl = 1/(1+mu_t);
            dUp_dysl = 1./(1./(mu_half/muw.*dU_incomp_dysl)-f);
            % above should use avg of [rho, mu, y_plus] since all the derivs are
            % centered and then midpoint rule is recovered
            u_sol(i_y) = u_sol(i_y-1) + dysl*dUp_dysl*utau;
            dT_du_w = s_Pr_input*(Tr-Tw)/Ue;
            T_next = Tw+dT_du_w*(u_sol(i_y)-u_sol(i_y).^2/u1)+(u_sol(i_y)/u1).^2.*(T1-Tw);
            T_sol(i_y) = T_prev + relax_T*(T_next-T_prev);
            T_sol(i_y) = max(T_sol(i_y),0);
            % check if converged
            if (T_sol(i_y)-T_prev)^2/T_sol(i_y)^2 < tol_T
                break;
            end
            T_prev = T_sol(i_y);
        end % end temperature iteration
        if iter_sec == max_iter_sec
            fprintf('\tT unconverged after max_iter_sec of %i reached. res=%f\n',iter_sec,(T_sol(i_y)-T_prev)^2/T_sol(i_y)^2)
        end
        if i_y+1 >= Ny
            assert(0); % need to allocate more memory
        end
        rho_m1 = rho;
        mu_m1 = mu;
        y_sl_m1 = y_sl;
        
        dy = max(dy,step_ctl_fax*u_sol(i_y)/abs(dysl*dUp_dysl*utau/dy));

    end % end while y < y1
    assert(abs(y_sol(i_y)-y1)/y1 < 1.0e-12)

    % secant method iteration:
    if iter == 1 % special treatment for first step
        utau_next = utau*1.1;
    else
        utau_next = utau + relax_u*(u1-u_sol(i_y))*(utau-utau_prev)/(u_sol(i_y)-U_end_prev);
    end
    if utau_next == inf
        utau_next = utau*0.9;
    end
    if utau_next < 0
        utau_next = utau*1.1;
    end
    U_end_prev = u_sol(i_y);
    utau_prev = utau;
    if abs(u1-u_sol(i_y))/u1 < tol_u
        fprintf('utau converged in %i iterations\n',iter)
        utau = utau_next;
        tau_w = rhow*utau^2;
        break;
    end
    utau = utau_next;
    tau_w = rhow*utau^2;
end % end nonlinear iteration of the equations
if iter == max_iter
    fprintf('utau unconverged after max_iter of %i reached. res=%f\n',iter,abs(utau_next-utau)/utau)
end
y_sol = y_sol(1:i_y);
u_sol = u_sol(1:i_y);
T_sol = T_sol(1:i_y);

% use the Reynolds analogy
qw = s_Pr_input/Pr*tau_w*Cp*(Tw-Tr)/Ue;
end % end of function
