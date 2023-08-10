clear all; 
close all;
set(0,'defaultAxesFontSize',18)
set(0,'defaulttextInterpreter','latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultTextInterpreter', 'latex')

% Select which case to run by uncommenting 1 of the following lines:
% case_1 = load('reference_data/cbl_dns_zhang_5.mat');
case_1 = load('reference_data/cbl_dns_volpiani_8.mat');
% case_1 = load('reference_data/cchannel_dns_trettel_4.mat');
% case_1 = load('reference_data/cchannel_dns_trettel_5.mat');
% case_1 = load('reference_data/cchannel_dns_trettel_8.mat');
% case_1 = load('reference_data/cchannel_dns_trettel_9.mat');
% case_1 = load('reference_data/cpipe_dns_modesti_3.mat');
% case_1 = load('reference_data/cpipe_dns_modesti_4.mat');

% Assign figure number to the plots.
plot_U = 1;
plot_T = 2;
figure(plot_U); hold on;
figure(plot_T); hold on;

% Load incompressible reference data.
incomp_data = load('reference_data/channel_dns_lee_5.mat');

%% Load the present i_case into the matlab workspace
cellfun(@(x,y) assignin('base',x,y),fieldnames(case_1),struct2cell(case_1));
Retau_star = y_sl(i_edge);
u_sl = sqrt(tau_w./(bar_rho_rhow*rhow));
l_sl = bar_mu_muw*muw./(bar_rho_rhow*rhow)./u_sl;
y_dim = y_plus*l_visc;
fprintf('%s: Retau=%.1f, Retau_star=%.1f, Me=%.4f, Minf=%.4f, Mb=%.2f\n',case_type ,Retau, Retau_star, Me, bar_M(end), Mb)
fprintf('Re_b=%.2f, Re_bw=%.2f\n',rhob*Ub*del99/muw,rhow*Ub*del99/muw)
Cp = gamma*R/(gamma -1); % perfect gas
qw = -Cp*muw/Pr*(bar_T_Tw(2)-bar_T_Tw(1))*Tw/(y_dim(2)-y_dim(1));
fprintf('B_q=%.3f\n',qw/(rhow*Cp*utau*Tw))

% Plot compressible DNS reference data.
figure(plot_U); plot(y_plus,bar_U_utau,'-k','displayname','DNS','linewidth',1.5)
figure(plot_T); plot(y_del99,bar_T_Tw,'-k','displayname','DNS','linewidth',1.5)

% Choose the matching location.
y1 = 0.3*del99; % the matching location in dimensional units
u1 = interp1(y_del99,U_utau*utau,y1/del99,'linear','extrap');
rho1 = interp1(y_del99,bar_rho_rhow*rhow,y1/del99,'linear','extrap');
mu1 = interp1(y_del99,bar_mu_muw*muw,y1/del99,'linear','extrap');
T1 = interp1(y_del99,bar_T_Tw*Tw,y1/del99,'linear','extrap');
Cp = gamma*R/(gamma -1); % perfect gas
Pr_t = 0.9;
r = 0.89;
% Set up constant for viscosity law.
if strcmp(visc_law,'power law')
    visc_arg = visc_power;
elseif strcmp(visc_law,'sutherland')
    visc_arg = S_sutherland_Tw*Tw;
else
    assert(0);
end

Tr_GRA = Te*(1+r*(gamma-1)/2*Me^2);
qw = -Cp*muw/Pr*(bar_T_Tw(2)-bar_T_Tw(1))*Tw/(y_dim(2)-y_dim(1));
s_Pr_DNS = qw*Ue/(tau_w*Cp*(Tw-Tr_GRA))*Pr;

%% Classical ODE model
[u_sol_co, y_sol_co, tau_w_co, T_sol_co, qw_co] = classical_wm(y1,u1,rho1,mu1,T1,Tw,case_1,Pr,Pr_t,visc_law,visc_arg);
utau_co = sqrt(tau_w_co/rhow);
u_plus_co = u_sol_co/utau_co;
y_plus_co = y_sol_co/muw*utau_co*rhow;

%% Compute wall modeled profiles
s_DM4 = 1.14;
r = Pr^(1/3);
s_Pr = s_DM4*Pr;
[u_sol_shoot, y_sol_shoot, tau_w_shoot, T_sol_shoot, qw_shoot] = inv_transf_wm(y1,u1,rho1,mu1,T1,Tw,case_1,r,visc_law,visc_arg,s_Pr);
utau_shoot = sqrt(tau_w_shoot/rhow);
u_plus_shoot = u_sol_shoot/utau_shoot;
y_plus_shoot = y_sol_shoot/muw*utau_shoot*rhow;

figure(plot_U);
plot(y_plus_co,u_plus_co,'--b','displayname','Classical WM','linewidth',2)
plot(y_plus_shoot,u_plus_shoot,'-.r','displayname','Present WM','linewidth',2)

figure(plot_T); 
plot(y_sol_co/del99,T_sol_co/Tw,'--b','displayname','Classical WM','linewidth',2)
plot(y_sol_shoot/del99,T_sol_shoot/Tw,'-.r','displayname','Present WM','linewidth',2)

%% Formatting the plots
figure(plot_U)
set(gca,'xscale','log')
xlim([1 Retau])
ylab = ylabel('$\tilde{U}^+$');
xlab = xlabel('$y^+$');
legend('location','best')
legend('boxoff'); box on;

figure(plot_T)
set(gca,'xscale','log')
xlim([1/Retau 1])
ylab = ylabel('$\tilde{T}/\tilde{T}_w$');
xlab = xlabel('$y/\delta$');
legend('location','best')
legend('boxoff'); box on;
