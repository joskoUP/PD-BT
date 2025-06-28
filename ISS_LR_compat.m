% This is a script to create plots for PD-BT for the ISS1R model with a
% compatible prior covariance.
% Some of the code for balancing Bayesian inference / Spantini is taken
% from https://github.com/elizqian/balancing-bayesian-inference.
% (Requires lyapchol from control system toolbox).

clear; close all
rng(1)
%% define LTI model
% ISS model
load('iss1R.mat')
d           = size(A,1);
d_out       = size(C,1);
A_full      = A;
C_full      = C;
sig_obs     = [2.5e-3 5e-4 5e-4]';

%% given full low-rank prior - sample covariance of compatible covariance
ensemble_size       = 90; % size of prior ensemble
Lyap_solution       = lyapchol(A_full,B)'; % compatibility Lyapunov eq.
ensemble            = Lyap_solution * randn(d,ensemble_size); 
% ensemble drawn from compatible covariance 
Gamma_pr_full       = cov(ensemble'); % prior is ensemble covariance

% check and force prior compatibility
L_pr_full       = minmodPriorCompat(Gamma_pr_full,A_full)';
Gamma_pr_full   = L_pr_full*L_pr_full';

prior_rank      = rank(Gamma_pr_full);


% compatibility check
compat      = A_full*Gamma_pr_full+Gamma_pr_full*A_full';
eigcompat   = eig(compat);
eigcompat   = eigcompat/max(abs(eigcompat));
if sum(real(eigcompat)>1e-12) > 0
    disp('The prior covariance L_pr is not prior compatible.')
end

%% compute full noisy observability Gramian
% helper matrix    
F           = C_full./sig_obs;
L_Q         = lyapchol(A_full',F')';
Q_inf       = L_Q*L_Q';

%% compute prior-driven reachability Gramian
R_P         = lyapchol(A_full,L_pr_full)';  
P           = R_P*R_P';

%% set up the experiment
% generate random data from multiple initial conditions
num_reps    = 100;
x0_all      = L_pr_full*randn(size(L_pr_full,2),num_reps);

% define time frame (measurement times) for the inference problem
T               = 8; % end time
dt_obs          = 0.1; % difference between measurements   
n               = round(T/dt_obs); % number of measurements
sig_obs_long    = repmat(sig_obs,n,1);

% define full forward model 
G       = zeros(n*d_out,d);
iter    = expm(A_full*dt_obs);
temp    = C;
for i   = 1:n
    temp                        = temp*iter;
    G((i-1)*d_out+1:i*d_out,:)  = temp;
end

% compute Fisher info
Go  = G./sig_obs_long;
H   = Go'*Go;
% helper matrix
M   = G*L_pr_full;

% generate measurements - multiple measurements for Bayes risk
y_all   = G*x0_all;
m_all   = y_all + sig_obs_long.*randn(n*d_out,num_reps);

%% compute true posterior
full_rhs_all        = G'*(m_all./(sig_obs_long.^2));
PosCov_true         = L_pr_full*(eye(size(L_pr_full,2)) - M'*((M*M'+diag(sig_obs_long.^2))\M))*L_pr_full';
PosMean_true_all    = PosCov_true*full_rhs_all;

%% compute restricted full-rank model for space where prior is invertible
% we use the prior_rank most important directions of Gamma_pr*H
[U,S,Z]             = svd(Go*L_pr_full); 
Gamma_pr            = S(1:prior_rank,1:prior_rank);
tau_full            = diag(S);
L_pr                = sqrt(Gamma_pr);
L_prinv             = diag(1./diag(L_pr));
Vg                  = (L_prinv*U(:,1:prior_rank)'*Go)';
Wg                  = L_pr_full*Z(:,1:prior_rank)*L_prinv;
% non-normalized balancing transformation for restriction
Wtilde_full         = (L_prinv*L_prinv*U(:,1:prior_rank)'*Go)';
What_full           = L_pr_full*Z(:,1:prior_rank);
% normalized Spantini directions

% restricted forward matrix
G_res               = G*Wg; % restricted forward matrix
Go_res              = G_res./sig_obs_long; % restricted sqrt factor of H

% restricted posterior quantities (we now can invert the restricted prior)
R_posinv            = qr([Go_res; L_prinv],0);
R_posinv            = triu(R_posinv(1:prior_rank,:)); % Pull out upper triangular factor
L_pos_res           = inv(R_posinv);
PosCov_res          = L_pos_res*L_pos_res';
full_rhs_res_all    = G_res'*(m_all./(sig_obs_long.^2));
PosMean_res_all     = R_posinv\(R_posinv'\full_rhs_res_all);

%% compute posterior approximations and errors
r_vals      = 1:5:prior_rank; % reduced rank must be smaller than prior_rank
rmax        = max(r_vals);

% computations with restricted (H,Gamma_pr^-1) (Spantini LIS)
[~,S,W]     = svd(Go_res*L_pr);
tau         = diag(S);
What        = L_pr*W;       % spantini directions
Wtilde      = L_prinv*W;

% balancing with (Q_m,Gamma_pr^-1) (LIS_BT) in full space
[V,S,W]     = svd(L_Q'*L_pr_full); 
S           = S(1:rmax,1:rmax);
delQ        = diag(S);
Siginvsqrt  = diag(1./sqrt(delQ));
Sr          = (Siginvsqrt*V(:,1:rmax)'*L_Q')';
Tr          = L_pr_full*W(:,1:rmax)*Siginvsqrt; % balancing transformation
A_LIS_BT    = Sr'*A_full*Tr;
C_LIS_BT    = C_full*Tr;

% prior-driven balancing (PD) with (Q_m,P) with A_fullP=PA_full'=-Gamma_pr 
[V,S,W]     = svd(L_Q'*R_P); 
% compute balanced quantities for posterior computation and error bound
delQ        = diag(S);
Siginvsqrt  = diag(1./sqrt(delQ));
Sr_PD       = (Siginvsqrt*V'*L_Q')';
Tr_PD       = R_P*W*Siginvsqrt; % balancing transformation
A_PD        = Sr_PD'*A_full*Tr_PD;
C_PD        = C_full*Tr_PD;
L_pr_PD     = Sr_PD'*L_pr_full;

%% compute posterior approximations
f_dist          = zeros(length(r_vals),3); % Förstner metric (restricted space)
F_dist          = zeros(length(r_vals),3); % relative Frobenius error (full space)
mu_errs         = zeros(length(r_vals),3); % empirical Bayes risk (restricted space)
full_mu_errs    = zeros(length(r_vals),3); % mean square error (full space)
error           = zeros(length(r_vals),1); % output error full system vs. PD-BT
Trace_vals      = zeros(length(r_vals),1); % trace for PD-BT error bound
kappa           = zeros(length(r_vals),1); % kappa for PD-BT error bound (h)
for rr = 1:length(r_vals)
    r           = r_vals(rr);
        
    %% Spantini posterior quantities
    % Spantini approx posterior covariance in restricted space
    PosCov_sp       = What*diag([1./(1+tau(1:r).^2); ones(prior_rank-r,1)])*What'; 
    % Spantini posterior in full space
    PosCov_sp_full  = What_full*diag([1./(1+tau_full(1:r).^2); ones(prior_rank-r,1)])*What_full'; 
    f_dist(rr,1)    = forstner(PosCov_sp,PosCov_res,'spd');
    F_dist(rr,1)    = norm(PosCov_sp_full-PosCov_true,'fro')/norm(PosCov_true,'fro');
    
    % Spantini approx posterior mean in restricted space
    Pi_r                = What(:,1:r)*Wtilde(:,1:r)'; % Spantini projector
    mean_sp_1           = PosCov_sp*Pi_r'*full_rhs_res_all;
    temp_sp_1           = L_pos_res\(mean_sp_1 - PosMean_res_all);
    mu_errs(rr,1)       = mean(sqrt(sum(temp_sp_1.^2)));
    Pi_r_full           = What_full(:,1:r)*Wtilde_full(:,1:r)'; % Spantini projector
    mean_sp_1           = PosCov_sp_full*Pi_r_full'*full_rhs_all;
    full_mu_errs(rr,1)  = mean(sum((mean_sp_1-PosMean_true_all).^2));
    % mean square error for mean approximation in full space
    
    %% LIS_BT posterior quantities
    % LIS_BT - generate reduced forward matrix G_LIS_BT (not restricted)
    G_LIS_BT        = zeros(n*d_out,r);
    iter            = expm(A_LIS_BT(1:r,1:r)*dt_obs);
    temp            = C_LIS_BT(:,1:r);
    for i = 1:n
        temp                            = temp*iter;
        G_LIS_BT((i-1)*d_out+1:i*d_out,:)  = temp;
    end
    G_LIS_BT        = G_LIS_BT*Sr(:,1:r)';
    M_LIS_BT        = G_LIS_BT*L_pr_full; % helper matrix for posterior
   
    % LIS_BT - compute posterior covariance
    PosCov_LIS_BT       = L_pr_full*(eye(size(L_pr_full,2)) - M_LIS_BT'*((M_LIS_BT*M_LIS_BT'+diag(sig_obs_long.^2))\M_LIS_BT))*L_pr_full';
    PosCov_LIS_BT_res   = Vg'*PosCov_LIS_BT*Vg; % restricted LIS_BT posterior covariance
    f_dist(rr,2)        = forstner(PosCov_LIS_BT_res,PosCov_res,'spd');
    F_dist(rr,2)        = norm(PosCov_LIS_BT-PosCov_true,'fro')/norm(PosCov_true,'fro');
    % LIS_BT - compute posterior mean
    mean_LIS_BT         = PosCov_LIS_BT*G_LIS_BT'*(m_all./(sig_obs_long.^2));
    temp_LIS_BT         = L_pos_res\(Vg'*mean_LIS_BT - PosMean_res_all);
    mu_errs(rr,2)       = mean(sqrt(sum(temp_LIS_BT.^2))); % Bayes risk error in restricted space
    full_mu_errs(rr,2)  = mean(sum((mean_LIS_BT-PosMean_true_all).^2));

    %% PD-BT posterior quantities
    % PD-BT - generate reduced forward matrix G_LIS_BT (not restricted)
    G_PD           = zeros(n*d_out,r);
    iter           = expm(A_PD(1:r,1:r)*dt_obs);
    temp           = C_PD(:,1:r);
    for i = 1:n
        temp                             = temp*iter;
        G_PD((i-1)*d_out+1:i*d_out,:)   = temp;
    end
    G_PDS          = G_PD*Sr_PD(:,1:r)';
    M_PD           = G_PD*L_pr_PD(1:r,:); % reduced helper matrix
    
    % PD-BT - compute posterior covariance
    PosCov_PD      = L_pr_full*(eye(size(L_pr_full,2)) - M_PD'*((M_PD*M_PD'+diag(sig_obs_long.^2))\M_PD))*L_pr_full';
    PosCov_PD_res  = Vg'*PosCov_PD*Vg; % restricted PD-BT posterior covariance    
    f_dist(rr,3)   = forstner(PosCov_PD_res,PosCov_res,'spd');
    F_dist(rr,3)   = norm(PosCov_PD-PosCov_true,'fro')/norm(PosCov_true,'fro');
    % PD-BT - compute posterior mean
    mean_PD        = PosCov_PD*G_PDS'*(m_all./(sig_obs_long.^2));
    temp_PD        = L_pos_res\(Vg'*mean_PD - PosMean_res_all);
    mu_errs(rr,3)  = mean(sqrt(sum(temp_PD.^2))); % Bayes risk error in restricted space
    full_mu_errs(rr,3)   = mean(sum((mean_PD-PosMean_true_all).^2));

    % PD-BT - assemble necessary quantities for the outpur error bound
    y_all_PD       = G_PDS*x0_all;
    error(rr)      = mean(sum(((y_all-y_all_PD)./sig_obs_long).^2));
    L2             = L_pr_PD(r+1:d,:);
    C1             = C_PD(:,1:r)./sig_obs;
    A11            = A_PD(1:r,1:r);
    A12            = A_PD(1:r,r+1:d);
    Y              = sylvester(A_PD',A11,-C_PD'*C1);
    Y2             = Y(r+1:d,:);
    delQ2          = diag(delQ(r+1:d));
    Trace_vals(rr) = trace((L2*L2'+2*Y2*A12)*delQ2);

    % compute kappa for the  error bound based on the impulse response
    % this takes very long compared to the rest of the script
    % for quick PD-BT runs, comment out and use kappa = 10
    fun = @(t) norm(((F*expm(A*t)*L_pr_full)-((C_PD(:,1:r)./sig_obs)*expm(A_PD(1:r,1:r)*t)*L_pr_PD(1:r,:))),'fro')^2;
    C = zeros(1,n);
    for i=1:n
        t = i*dt_obs;
        left    = fun(t);
        right   = integral(fun,t-dt_obs,t,'ArrayValued',true);
        C(i)    = left/right;
    end
    kappa(rr) = max(C);
end

%% plots
% Warning if complex parts of Förstner distances are nontrivial
if ~isempty(find(abs(imag(f_dist))>eps*abs(real(f_dist)), 1))
    warning('Significant imaginary parts found in Forstner distance')
end
% Otherwise imaginary parts are trivial artifacts of generalized eig
f_dist  = real(f_dist);

figure; clf  
%% plot posterior covariance Forstner errors (in restricted space)
subplot(1,2,2)
semilogy(r_vals,f_dist(:,2),'o','Color','#FFB000'); hold on
semilogy(r_vals,f_dist(:,3),'^','Color','#DC267F');
semilogy(r_vals,f_dist(:,1),'x-','Color','#648FFF','LineWidth',2);
grid on
title('$\hat{\mathbf{\Gamma}}_{pos}$ -- F\"orstner error','interpreter','latex','fontsize',20) 
xlabel('$r$','interpreter','latex','fontsize',13)
legend({'LIS-BT','PD-BT','OLR'},'interpreter','latex','fontsize',13,'Location','northeast')
legend boxoff
xlim([0 236])
ylim([1e-12 1e+4])
set(gca,'fontsize',13,'ticklabelinterpreter','latex')

%% plot posterior mean errors in relative PosCov^-1 norm (in restricted space)
posnormref = mean(sqrt(sum((L_pos_res\(PosMean_res_all)).^2)));
subplot(1,2,1)
semilogy(r_vals,mu_errs(:,2)/posnormref,'o','Color','#FFB000'); hold on
semilogy(r_vals,mu_errs(:,3)/posnormref,'^','Color','#DC267F');
semilogy(r_vals,mu_errs(:,1)/posnormref,'x-','Color','#648FFF','LineWidth',2);
grid on
title('$\hat{\mathbf{\mu}}_{pos}$ -- Bayes risk','interpreter','latex','fontsize',20)
xlabel('$r$','interpreter','latex','fontsize',13)
legend({'LIS-BT','PD-BT','OLR'},'interpreter','latex','fontsize',13,'Location','northeast')
legend boxoff
xlim([0 236])
ylim([1e-12 1])
set(gca,'fontsize',13,'ticklabelinterpreter','latex')
   
figure; clf 
%% plot posterior covariance Frobenius errors (in full space)
subplot(1,2,2)
semilogy(r_vals,F_dist(:,2),'o','Color','#FFB000'); hold on
semilogy(r_vals,F_dist(:,3),'^','Color','#DC267F');
semilogy(r_vals,F_dist(:,1),'x-','Color','#648FFF','LineWidth',2);
grid on
title('$\mathbf{\Gamma}_\mathrm{pos}$ -- relative Frobenius error','interpreter','latex','fontsize',20) 
xlabel('$r$','interpreter','latex','fontsize',13)
legend({'LIS-BT','PD-BT','OLR'},'interpreter','latex','fontsize',13,'Location','northeast')
legend boxoff
xlim([0 236])
ylim([1e-11 1e+3])
set(gca,'fontsize',13,'ticklabelinterpreter','latex')

%% plot posterior mean square errors (in full space)
posnormref = mean(sqrt(sum((PosMean_true_all).^2)));
subplot(1,2,1)
semilogy(r_vals,full_mu_errs(:,2)/posnormref,'o','Color','#FFB000'); hold on
semilogy(r_vals,full_mu_errs(:,3)/posnormref,'^','Color','#DC267F');
semilogy(r_vals,full_mu_errs(:,1)/posnormref,'x-','Color','#648FFF','LineWidth',2);
grid on
title('$\mathbf{\mu}_\mathrm{pos}$ -- mean square error','interpreter','latex','fontsize',20)
xlabel('$r$','interpreter','latex','fontsize',13)
legend({'LIS-BT','PD-BT','OLR'},'interpreter','latex','fontsize',13,'Location','northeast')
legend boxoff
xlim([0 236])
ylim([1e-12 100])
set(gca,'fontsize',13,'ticklabelinterpreter','latex')

%% plot PD-BT error bound    
figure; clf
semilogy(r_vals,error,'o-','Color','#DC267F','LineWidth',2); hold on
semilogy(r_vals,prior_rank*kappa.*Trace_vals,'x-','Color','#648FFF','LineWidth',2)
title('Comparison between output error and bound','interpreter','latex','fontsize',20)
xlabel('$r$','interpreter','latex','fontsize',13)
legend({'$E_{\mathbf{p}}\|\mathbf{Y}-\hat{\mathbf{Y}}\|_{\mathbf{\Gamma}_\epsilon^{-1}}^2$','error bound'},'interpreter','latex','fontsize',13,'Location','southwest')
legend boxoff
xlim([0 prior_rank])
ylim([1e-15 1e+7])
set(gca,'fontsize',13,'ticklabelinterpreter','latex')
