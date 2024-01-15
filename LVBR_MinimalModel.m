% LV Corticothalamic Model of Visual Rivalry

% Christopher Whyte 

rng('shuffle')
close all; clear

%% Parameters  

%----- simulation length

sim_time = 20000;

%----- number of neurons 

n_e = 90;     % excitatory LV cortical neurons 
n_i = n_e;    % inhibitory LV cortical interneurons 
n_th = 10;    % non-specific matrix thalamus neurons

%----- connectivity 

% E/I ring 
e_e = 3.25; e_i = 1; i_e = 4.5; 
sd_e_e = .5; sd_e_i = 2; sd_i_e = 2; 
J_e_e = zeros(n_e,n_e); J_e_i = zeros(n_e,n_e); J_i_e = zeros(n_e,n_e); 

% Cortical ring coupling
theta = linspace(0,2*pi,n_e);  % location on ring
for j = 1:n_e
    for k = 1:n_e
        % euclidean distance in coordinates of unit circle
        dist = sqrt(((cos(theta(j)) - cos(theta(k)))^2 + (sin(theta(j)) - sin(theta(k)))^2)); 
        J_e_e(j,k) = e_e*exp(-.5*(dist./sd_e_e).^2)./(sd_e_e*sqrt(2*pi)); 
        J_e_i(j,k) = e_i*exp(-.5*(dist./sd_e_i).^2)./(sd_e_i*sqrt(2*pi));
        J_i_e(j,k) = i_e*exp(-.5*(dist./sd_i_e).^2)./(sd_i_e*sqrt(2*pi));
    end 
end


% ring -> th
e_th = 4;
J_e_th = zeros(n_th,n_e);
pocket_size = 9; counter = 1; 
for j = 1:n_th
    J_e_th(j,counter:counter+pocket_size-1) = e_th;
    counter = counter + pocket_size;
end

% ring <- th 
th_d = 1;
J_th_d = zeros(n_th,n_e);
pocket_size = 9; counter = 1; 
for j = 1:n_th
    J_th_d(j,counter:counter+pocket_size-1) = th_d;
    counter = counter + pocket_size;
end

%----- neuron params 

% C = capicitance
% v_r = rest potential
% v_th = threshold
% v_peak = spike cut off
% k = membrane potential time constant
% a = adaptation time constant
% b = sensitivity to subthreshold oscillations
% c = voltage reset
% d = adaptation reset 

% L5b dendrite params 
vd_peak = 20; 
Cd = 100; vd_r = -50; vd_t = -50; kd = 1;
ad = .8; bd = 15; cd = -10; dd = 500;

% L5b soma params 
ve_peak = 50; 
Ce = 150; ve_r = -75; ve_t = -45; ke = 1.2;
ae = 0.01; be = 5; 
ce = -66*ones(n_e,sim_time); burst_c_reset = -56;
de = 260*ones(n_e,sim_time); burst_d_reset = 130;

% basket cell params
vi_peak = 25; 
Ci = 20; vi_r = -55; vi_t = -40; ki = 1;
ai = 0.15; bi = 8; ci = -55; di = 200;

% matrix thalamus params 
vth_peak = 35;
C_th = 200; vth_r = -60; vth_t = -50; kth = 1.6; 
ath = 0.01; bth = 15; cth = -60; dth = 10; 

%----- synapse params

% reverse potentials
E_exc = 0; E_inh = -75; E_adapt = -80; 

% adaptation current
tau_decay_adapt = 2000; delta_adapt = 0.045;

% NMDA
beta = 1/(3.57); alpha = 0.062; Mg = 1.2;
tau_rise_NMDA = 10;  tau_decay_NMDA = 75;
NMDA_window = 3000;

% AMPA
tau_rise_AMPA = 1;   tau_decay_AMPA = 8; 
AMPA_window = 500;

% GABAa
tau_rise_GABAa = 1;  tau_decay_GABAa = 6; 
GABAa_window = 500;

%----- input into ring

% poisson process drive to soma
fr_L = 475;                           % L firing rate in hertz
fr_R = 475;                           % R firing rate in hertz
tsim = sim_time/1000;                 % simulation length (seconds)
dt = 1/1000;                          % time step
nbins = floor(tsim/dt);               % time bins
spikemat_L = (rand(1,nbins)<fr_L*dt); % left eye spikes
spikemat_R = (rand(1,nbins)<fr_R*dt); % right eye spikes

% FF synaptic input (convolution of exponential synapses with spike train)
synapse = @(syn_window,tau_rise,tau_decay,spikes) conv(spikes,(1-exp(-(1:syn_window)./tau_rise)).*exp(-(1:syn_window)./tau_decay));

% left eye
g_syn_AMPA_L = synapse(sim_time, tau_rise_AMPA, tau_decay_AMPA,spikemat_L); g_syn_AMPA_L = g_syn_AMPA_L(1:sim_time);
g_syn_NMDA_L = synapse(sim_time, tau_rise_NMDA, tau_decay_NMDA,spikemat_L); g_syn_NMDA_L = g_syn_NMDA_L(1:sim_time);

% right eye
g_syn_AMPA_R = synapse(sim_time, tau_rise_AMPA, tau_decay_AMPA,spikemat_R); g_syn_AMPA_R = g_syn_AMPA_R(1:sim_time);
g_syn_NMDA_R = synapse(sim_time, tau_rise_NMDA, tau_decay_NMDA,spikemat_R); g_syn_NMDA_R = g_syn_NMDA_R(1:sim_time);

NL = 23; NR = NL+.5*n_e; sd = 10; 
I_ext_const = zeros(n_e,1); 
I_ext_NMDA = zeros(n_e,sim_time); 
I_ext_AMPA= zeros(n_e,sim_time);
for n = 1:n_e
    I_ext_const(n,1) = 200.*(exp(-((n-NL)/sd).^2) + exp(-((n-NR)/sd).^2));
    I_ext_NMDA(n,:) = exp(-((n-NL)/sd).^2).*g_syn_NMDA_L + exp(-((n-NR)/sd).^2).*g_syn_NMDA_R;
    I_ext_AMPA(n,:) = exp(-((n-NL)/sd).^2).*g_syn_AMPA_L + exp(-((n-NR)/sd).^2).*g_syn_AMPA_R;
end 


%% Model Lesions/Stimulation

Baclofen = zeros(n_e,1);
% Baclofen(1:45,:) = Baclofen(1:45,:) + 100;
% Baclofen(1:45,:) = Baclofen(1:45,:) - 100;

%% Initialise 

%----- Initialise state variable matrices
vd = zeros(n_e,sim_time); ud = zeros(n_e,sim_time);
ve = zeros(n_e,sim_time); ue = zeros(n_e,sim_time);
vi = zeros(n_i,sim_time); ui = zeros(n_i,sim_time);
vth = zeros(n_th,sim_time); uth = zeros(n_th,sim_time);

%----- Initialise synapse matrices

g_syn_e_e_NMDA = zeros(n_e,sim_time);
g_syn_e_e_AMPA = zeros(n_e,sim_time);
g_syn_e_i_AMPA = zeros(n_i,sim_time);
g_syn_e_i_NMDA = zeros(n_i,sim_time);
g_syn_e_th_AMPA = zeros(n_th,sim_time);
g_syn_i_e_GABAa = zeros(n_e,sim_time);
g_syn_th_d_NMDA = zeros(n_e,sim_time);
g_syn_th_d_AMPA = zeros(n_e,sim_time);
g_syn_adapt = zeros(n_e,sim_time);

%----- Initial conditions

vd(:,1) = vd_r;
ve(:,1) = ve_r; 
vi(:,1) = vi_r;
vth(:,1) = vth_r;

%----- initialise spike storage containers

firings_d = []; 
firings_e = []; tf_e = zeros(n_e,1); 
firings_i = []; tf_i = zeros(n_i,1); 
firings_th = []; tf_th = zeros(n_th,1); 

%% Simulation loop

disp('Simulation start'); tic
for t = 1:(sim_time-1)
    
    % display sim time
    if mod(t,1000)==0
        disp(['t = ',num2str(t)]);
    end
    
    %---------- L5b apical dendrites
    
    % conductance
    Mg_gate = (1./(1 + Mg.*beta.*exp(-alpha*vd(:,t))));
    g_syn_tot = Mg_gate.*g_syn_th_d_NMDA(:,t) + g_syn_th_d_AMPA(:,t) + eps; % total condunctance
    E_tot = (Mg_gate.*g_syn_th_d_NMDA(:,t).*E_exc + g_syn_th_d_AMPA(:,t).*E_exc )./g_syn_tot; % total reverse potential
    % membrane potential
    vd(:,t+1) = (vd(:,t) + (kd*(vd(:,t) - vd_r).*(vd(:,t) - vd_t) - ud(:,t) + 20*(ve(:,t) - vd(:,t)) + Baclofen + g_syn_tot.*E_tot)./Cd)...
                ./(1 + g_syn_tot./Cd);    
    fired=vd(:,t+1)>=vd_peak;
    if any(fired==1) 
       firings_d = [firings_d; t+0*find(fired==1), find(fired==1)]; 
       % interpolate spike time to find tau
       t_peak = t + (vd_peak - vd(fired,t))./(vd(fired,t+1) - vd(fired,t));
       tau = t_peak - t;
       % update adapatation variable with new tau
       ud(fired,t+1) = ud(fired,t) + tau.*ad.*(bd.*(vd(fired,t) - vd_r) - ud(fired,t));
       % update adaptation variable
       ud(fired,t+1) = ud(fired,t+1) + dd; 
       % pad spike amplitude
       vd(fired,t) = vd_peak; 
       vd(fired,t+1) = cd;
       % reset conditions on soma conditional on apical spike
       ce(fired,t:t+25) = burst_c_reset; de(fired,t:t+25) = burst_d_reset;
    end 
    ud(~fired,t+1) = ud(~fired,t)+ad*(bd*(vd(~fired,t)-vd_r)-ud(~fired,t));
    
    %---------- L5b soma
    
    % conductance
    Mg_gate = (1./(1 + Mg.*beta.*exp(-alpha*ve(:,t))));
    g_syn_tot = g_syn_e_e_AMPA(:,t) + Mg_gate.*min(g_syn_e_e_NMDA(:,t),130) + g_syn_i_e_GABAa(:,t) + Mg_gate.*I_ext_NMDA(:,t) + I_ext_AMPA(:,t) + g_syn_adapt(:,t) + eps; % total condunctance
    E_tot = (g_syn_e_e_AMPA(:,t).*E_exc + Mg_gate.*min(g_syn_e_e_NMDA(:,t),130).*E_exc + g_syn_i_e_GABAa(:,t).*E_inh + Mg_gate.*I_ext_NMDA(:,t).*E_exc + I_ext_AMPA(:,t).*E_exc + g_syn_adapt(:,t).*E_adapt)./g_syn_tot; % total reverse potential
    % membrane potential
    ve(:,t+1) = (ve(:,t) + (ke.*(ve(:,t) - ve_r).*(ve(:,t) - ve_t) - ue(:,t) + I_ext_const + g_syn_tot.*E_tot)./Ce)...
                ./(1 + (g_syn_tot./Ce));
    g_syn_adapt(:,t+1) = g_syn_adapt(:,t) + (-g_syn_adapt(:,t)./tau_decay_adapt);
    % find spikes
    fired=ve(:,t+1)>=ve_peak;
    firedind = find(ve(:,t+1)>=ve_peak); 
    if any(fired==1) 
       % store spikes
       firings_e = [firings_e; t+0*firedind, firedind]; 
       % interpolate spike time to find tau
       t_peak = t + (ve_peak - ve(fired,t))./(ve(fired,t+1) - ve(fired,t));
       tau = t_peak - t;
       % update adapatation variable with new tau
       ue(fired,t+1) = ue(fired,t) + tau.*ae.*(be.*(ve(fired,t)-ve_r) - ue(fired,t));
       % update adaptation variable
       ue(fired,t+1) = ue(fired,t+1) + de(fired,t);
       % pad spike amplitude
       ve(fired,t) = ve_peak; 
       ve(fired,t+1) = ce(fired,t);
       % spike times 
       tf_e(fired) = t;
       % adaptation current 
       g_syn_adapt(fired,t+1) = g_syn_adapt(fired,t+1) + delta_adapt;
       % NMDA conductance window
       if t > (sim_time - NMDA_window) 
           NMDA_idx = sim_time - t;
       else
           NMDA_idx = NMDA_window;
       end 
       tdash_NMDA=(t+1:(t+NMDA_idx));
       % L5b -> L5b: NMDA
       g_syn_e_e_NMDA(:,t+1:(t+NMDA_idx)) = g_syn_e_e_NMDA(:,t+1:(t+NMDA_idx)) + J_e_e*((1-exp(-(tdash_NMDA-tf_e)./tau_rise_NMDA)).*exp(-(tdash_NMDA-tf_e)./tau_decay_NMDA).*fired);
       % L5b -> basket: NMDA
       g_syn_e_i_NMDA(:,t+1:(t+NMDA_idx)) = g_syn_e_i_NMDA(:,t+1:(t+NMDA_idx)) + J_e_i*((1-exp(-(tdash_NMDA-tf_e)./tau_rise_NMDA)).*exp(-(tdash_NMDA-tf_e)./tau_decay_NMDA).*fired);
       % AMPA conductance window
       if t > (sim_time - AMPA_window) 
           AMPA_idx = sim_time - t;
       else
           AMPA_idx = AMPA_window;
       end 
       tdash_AMPA=(t+1:(t+AMPA_idx));
       % L5b -> L5b: AMPA
       g_syn_e_e_AMPA(:,t+1:(t+AMPA_idx)) = g_syn_e_e_AMPA(:,t+1:(t+AMPA_idx)) + J_e_e*((1-exp(-(tdash_AMPA-tf_e)./tau_rise_AMPA)).*exp(-(tdash_AMPA-tf_e)./tau_decay_AMPA).*fired);
       % L5b -> basket: AMPA
       g_syn_e_i_AMPA(:,t+1:(t+AMPA_idx)) = g_syn_e_i_AMPA(:,t+1:(t+AMPA_idx)) + J_e_i*((1-exp(-(tdash_AMPA-tf_e)./tau_rise_AMPA)).*exp(-(tdash_AMPA-tf_e)./tau_decay_AMPA).*fired);
       % L5b -> thal
       g_syn_e_th_AMPA(:,t+1:(t+AMPA_idx)) = g_syn_e_th_AMPA(:,t+1:(t+AMPA_idx)) + J_e_th*((1-exp(-(tdash_AMPA-tf_e)./tau_rise_AMPA)).*exp(-(tdash_AMPA-tf_e)./tau_decay_AMPA).*fired);
    end 
    ue(~fired,t+1) = ue(~fired,t) + ae.*(be.*(ve(~fired,t)-ve_r) - ue(~fired,t));
    
    %---------- basket cell
    
    % conductance
    Mg_gate = (1./(1 + Mg.*beta.*exp(-alpha*vi(:,t))));
    g_syn_tot = Mg_gate.*g_syn_e_i_NMDA(:,t) + g_syn_e_i_AMPA(:,t) + eps; % total condunctance
    E_tot = (Mg_gate.*g_syn_e_i_NMDA(:,t).*E_exc + g_syn_e_i_AMPA(:,t).*E_exc)./g_syn_tot; % total reverse potential
    % membrane potential
    vi(:,t+1) = (vi(:,t) + (ki.*(vi(:,t) - vi_r).*(vi(:,t) - vi_t) - ui(:,t) + g_syn_tot.*E_tot)./Ci)...
                ./(1 + (g_syn_tot./Ci));
    % find spikes
    fired=vi(:,t+1)>=vi_peak;
    firedind = find(vi(:,t+1)>=vi_peak); 
    if any(fired==1) 
       % store spikes
       firings_i = [firings_i; t+0*firedind, firedind]; 
       % interpolate spike time to find tau
       t_peak = t + (vi_peak - vi(fired,t))./(vi(fired,t+1) - vi(fired,t));
       tau = t_peak - t;
       % update adapatation variable with new tau
       ui(fired,t+1) = ui(fired,t) + tau.*ai.*(bi.*(vi(fired,t) - vi_r) - ui(fired,t));
       % update adaptation variable
       ui(fired,t+1) = ui(fired,t+1) + di;
       % pad spike amplitude
       vi(fired,t) = vi_peak; 
       vi(fired,t+1) = ci;
       % spike times 
       tf_i(fired) = t; 
       % GABAa conductance window
       if t > (sim_time - GABAa_window) 
           GABAa_idx = sim_time - t;
       else
           GABAa_idx = GABAa_window;
       end 
       tdash_GABAa=(t+1:(t+GABAa_idx));
       % basket -> L5b
       g_syn_i_e_GABAa(:,t+1:(t+GABAa_idx)) = g_syn_i_e_GABAa(:,t+1:(t+GABAa_idx)) + J_i_e*((1-exp(-(tdash_GABAa-tf_i)./tau_rise_GABAa)).*exp(-(tdash_GABAa-tf_i)./tau_decay_GABAa).*fired);
    end 
    ui(~fired,t+1) = ui(~fired,t) + ai.*(bi.*(vi(~fired,t) - vi_r) - ui(~fired,t));
    
    %---------- matrix thamalus population
    
    % conductance
    g_syn_tot = g_syn_e_th_AMPA(:,t) + eps; % total condunctance
    E_tot = (g_syn_e_th_AMPA(:,t).*E_exc)./g_syn_tot; % total reverse potential
    % membrane potential
    vth(:,t+1) = (vth(:,t) + (kth*(vth(:,t) - vth_r).*(vth(:,t) - vth_t) - uth(:,t) + g_syn_tot.*E_tot)./C_th)...
                 ./(1 + g_syn_tot./C_th);
    % find spikes
    fired=vth(:,t+1)>=vth_peak+0.1.*uth(t);
    if any(fired==1)
        firings_th = [firings_th; t+0*find(fired==1), find(fired==1)]; 
        % interpolate spike time to find tauth
        t_peak = t + (vth_peak - vth(fired,t))./(vth(fired,t+1) - vth(fired,t));
        tau = t_peak - t;
        % update adapatation vthariable with new tau
        uth(fired,t+1) = uth(fired,t) + tau.*ath.*(bth.*(vth(fired,t) - vth_r) - uth(fired,t));
        % update adaptation vthariable
        uth(fired,t+1) = uth(fired,t+1) + dth;
        % pad spike amplituthde
        vth(fired,t) = vth_peak; 
        % reset membrane potential
        vth(fired,t+1) = cth - 0.1.*uth(fired,t);
        % spike times 
        tf_th(fired) = t; 
        % NMDA conductance window
        if t > (sim_time - NMDA_window) 
           NMDA_idx = sim_time - t;
        else
           NMDA_idx = NMDA_window;
        end 
        tdash_NMDA=(t+1:(t+NMDA_idx));
        % thal -> dend NMDA
        g_syn_th_d_NMDA(:,t+1:(t+NMDA_idx)) = g_syn_th_d_NMDA(:,t+1:(t+NMDA_idx)) + J_th_d'*((1-exp(-(tdash_NMDA-tf_th)/tau_rise_NMDA)).*exp(-(tdash_NMDA-tf_th)/tau_decay_NMDA).*fired);
        % AMPA conductance window
        if t > (sim_time - AMPA_window) 
           AMPA_idx = sim_time - t;
        else
           AMPA_idx = AMPA_window;
        end 
        tdash_AMPA=(t+1:(t+AMPA_idx));
        % thal -> dend AMPA
        g_syn_th_d_AMPA(:,t+1:(t+AMPA_idx)) = g_syn_th_d_AMPA(:,t+1:(t+AMPA_idx)) + J_th_d'*((1-exp(-(tdash_AMPA-tf_th)/tau_rise_AMPA)).*exp(-(tdash_AMPA-tf_th)/tau_decay_AMPA).*fired);
    end 
    subthreshidx=vth(:,t+1) <= -65;
    if any(subthreshidx==1)
        uth(subthreshidx,t+1) = uth(subthreshidx,t) - ath*uth(subthreshidx,t);
    end 
    fired_subthresh = fired >= 1 | subthreshidx >= 1;
    uth(~fired_subthresh,t+1) = uth(~fired_subthresh,t) + ath*(bth*(vth(~fired_subthresh,t)-vth_r) - uth(~fired_subthresh,t));
    
end 
disp('Simulation end'); toc

%% basic figures

time = linspace(1,sim_time/1000,sim_time); % time in seconds
neuronidx = NL;
thalidx = round(neuronidx/9);

figure(1)
hold on
plot(time,ve(neuronidx,:),'k')
ylabel('Neuron index')
ylabel('Time (s)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')
title('L5 Soma')

figure(2)
plot(time,vi(neuronidx,:),'k') 
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')
title('Basket Cell')

figure(3)
plot(time,vth(thalidx,:),'k')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')
title('Thalamus')

figure(4)
hold on
title('Raster: L5 Soma')
plot(firings_e(:,1)/1000,firings_e(:,2),'.','Color', [220, 77 101]/255);
xlim([0,sim_time/1000])
ylim([0,n_e])
ylabel('Neuron index')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')

figure(5)
hold on
title('Raster: Basket Cell')
plot(firings_i(:,1)/1000,firings_i(:,2),'k.');
xlim([0,sim_time/1000])
ylim([0,n_i])
ylabel('Neuron index')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')

figure(6)
hold on
title('Raster: Thal')
plot(firings_th(:,1)/1000,firings_th(:,2),'k.');
xlim([0,sim_time/1000])
ylim([0,n_th])
ylabel('Neuron index')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')

figure(7)
plot(time,g_syn_adapt(neuronidx,:),'k')
title('Adaptation')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')

%% Average firing rate

firings = firings_e;

spikesLT = firings(1<= firings(:,2) & firings(:,2)<=45,1);
spikesRT = firings(45<=firings(:,2) & firings(:,2)<=90,1);

tstep = 1:1:sim_time;

rateL = histcounts(spikesLT,tstep)/45*1000; %count spikes in bin 1ms resol
rateR = histcounts(spikesRT,tstep)/45*1000;

SmoothGaussWindow = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');
rateL = SmoothGaussWindow(rateL,250);
rateR = SmoothGaussWindow(rateR,250);

figure(8)
hold on
plot(rateL,'k', 'linewidth', 2)
plot(rateR,'linewidth', 2,'Color', [220, 77 101]/255)
ylim([0,35])
ylabel('Firing rate (Hz)')
xlabel('Time (ms)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')




