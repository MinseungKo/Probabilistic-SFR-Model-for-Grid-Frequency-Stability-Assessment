clear all
close all
clc

%% Import Data
rng('shuffle')

PFile='loadbus_3_400.xlsx';
data = xlsread(PFile);
data = data(5:end, :);
time = data(:, 1);
Freq = data(:, 12:50);
Volt = data(:, 51:end);

PFile3='loadbus_7_500.xlsx';
data3 = xlsread(PFile3);
data3 = data3(5:end, :);
Freq3 = data3(:, 12:50);
Volt3 = data3(:, 51:end);

dist_bus_ref = 3;
dist_bus_tar = 7;
dP_ref = 400;
dP_tar = 500;

vsens = csvread('volt_sensitivity.csv'); %% Voltage sensitivity obtained by reactive power perturbation
H = csvread('voltH.csv');
rawFile = 'IEEE 39 bus_modified2.raw';

Ps = data(:, 2:11);
dt = time(2,1);
Tend = time(end,1);
N = round(Tend/dt)+1; % Total solve steps


%% Generator & System Informations
Gen_num = [2, 3, 4, 5, 6, 7, 8, 9, 10];
SG_num  = [2, 3, 4, 7, 9, 10];
Hs      = [2.3, 4.9, 4.3, 5.2, 3.15, 5.0];

Basepow = 100;
Sbase = sum(Ps(1,:))*Basepow;
Genbase = sum(Ps(1,SG_num))*Basepow;
GFLbase = Sbase-Genbase;

freqsys  = sum(Hs.*Freq(:,  SG_num),2)/sum(Hs);
freqsys3 = sum(Hs.*Freq3(:, SG_num),2)/sum(Hs);

Hsys = 0;
for i=1:6
    Hsys = Hsys+Hs(i)*Ps(1,SG_num(i))*100;
end
Dsys=0;
Hsys = Hsys/Sbase;

%% Generator Bus Voltage Estimation for Target Case
ini_load_v = [1.0345, 1.0221, 1.0102, 1.0098, 1.0109, 1.024, 1.0374, ...
                1.0355, 0.9941, 1.0359, 1.0474,1.0434,1.0586,1.0547,1.0428,1.0515,1.0518,0.982,1.03];
ini_load_p = [322, 400, 233.8, 522, 7.5, 220, 329.4, 158, 528, 274, 274.5, 208.6, 224, 139, 181, ...
                206,183.5, 9.2, 1104];

gen_bus = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39];
gen_load_bus = [2, 10];
load_bus = [3,4,7,8,12,15,16,18,20,21,23,24,25,26,27,28,29];

Vg_ref = Volt(:, gen_bus);
V_tar  = Volt3(:, gen_bus);
[T, Ng] = size(Vg_ref);

% pre-event samples
N0 = min(50, T);

V0_ref = mean(Vg_ref(1:N0,:), 1);
dV_ref = Vg_ref - V0_ref;

V0_tar = mean(V_tar(1:N0,:), 1);
dV_tar_true = V_tar - V0_tar;

t_event = time(N0);
t1 = t_event + 0.0;
t2 = 20.1;

% split point separating fast transient vs slow recovery
t_split = 8.0;        % tune: 6~10 sec
tau_bl  = 0.6;        % blend smoothness (sec)

idx_fast = find(time >= t1 & time <= t_split);
idx_slow = find(time >  t_split & time <= t2);

if numel(idx_fast) < 30
    idx_fast = (N0+1):min(T, N0+200);
end
if numel(idx_slow) < 30
    idx_slow = max(N0+1, round(0.5*T)):T;
end

% smooth blend gate 0->1
blend = 1./(1 + exp(-(time - t_split)/tau_bl));
blend(1:N0) = 0;
BlendT = repmat(blend, 1, Ng);

r_fast = 3;      % tune: 2~4
r_slow = 2;      % tune: 1~3

% FAST
Xf = dV_ref(idx_fast,:);
r_fast = min(r_fast, min(size(Xf)));
[Uf,Sf,Vf] = svd(Xf,'econ');
P_ref_fast = Vf(:,1:r_fast);
A_ref_fast = dV_ref * P_ref_fast / (P_ref_fast.'*P_ref_fast + 1e-12*eye(r_fast));

% SLOW
Xs = dV_ref(idx_slow,:);
r_slow = min(r_slow, min(size(Xs)));
[Us,Ss,Vs] = svd(Xs,'econ');
P_ref_slow = Vs(:,1:r_slow);
A_ref_slow = dV_ref * P_ref_slow / (P_ref_slow.'*P_ref_slow + 1e-12*eye(r_slow));


[dist_matrix, net] = build_dist_matrix_from_psse_raw(rawFile, 'X');

Nb = size(dist_matrix,1);
gen_buses_global = [30; 31; 32; 33; 34; 35; 36; 37; 38; 39];

% reduced Zbus to compute transfer factor g
Ybus = zeros(Nb,Nb);
for e = 1:numel(net.from)
    i = net.from(e); j = net.to(e);
    r = net.r(e); x = net.x(e);
    z = r + 1j*x;
    y = 1/z;
    Ybus(i,i) = Ybus(i,i) + y;
    Ybus(j,j) = Ybus(j,j) + y;
    Ybus(i,j) = Ybus(i,j) - y;
    Ybus(j,i) = Ybus(j,i) - y;
end

slack = 1;
idx = setdiff(1:Nb, slack);
Yred = Ybus(idx,idx);
Zred = inv(Yred);
bus2red = @(b) find(idx==b);

Zg_ref = zeros(Ng,1);
Zg_tar = zeros(Ng,1);
for ii = 1:Ng
    gi = gen_buses_global(ii);
    Zg_ref(ii) = abs( Zred(bus2red(gi), bus2red(dist_bus_ref)) );
    Zg_tar(ii) = abs( Zred(bus2red(gi), bus2red(dist_bus_tar)) );
end

g = (Zg_tar ./ (Zg_ref + 1e-12)) * (dP_tar/dP_ref);

% scalar location distance (for uncertainty scaling)
d_loc = dist_matrix(dist_bus_ref, dist_bus_tar);
if ~isfinite(d_loc) || d_loc<=0
    d_loc = 1.0;
end

Wg = zeros(Ng,Ng);
for i = 1:Ng
    for j = 1:Ng
        if i==j, continue; end
        Zi = bus2red(gen_buses_global(i));
        Zj = bus2red(gen_buses_global(j));
        zij = abs(Zred(Zi,Zj));
        Wg(i,j) = 1 / (zij + 1e-6);
    end
end
Wg = (Wg + Wg.')/2; Wg(1:Ng+1:end) = 0;
Dg = diag(sum(Wg,2));
Lg = Dg - Wg;

% diffusion operator for deterministic transfer
alpha_mix = 0.15;
Mmix = (eye(Ng) + alpha_mix*Lg) \ eye(Ng);

% gamma for fast modes
gamma_fast = zeros(r_fast,1);
for k = 1:r_fast
    pk = abs(P_ref_fast(:,k));
    pk = pk / (max(pk) + 1e-12);
    gamma_fast(k) = max(1 + max(pk), 0.8);
end

% gamma for slow modes
gamma_slow = zeros(r_slow,1);
for k = 1:r_slow
    pk = abs(P_ref_slow(:,k));
    pk = pk / (max(pk) + 1e-12);
    gamma_slow(k) = max(1 + max(pk), 0.8);
end

% FAST transfer
P_tar_fast = zeros(Ng, r_fast);
for k = 1:r_fast
    P_tar_fast(:,k) = Mmix * (P_ref_fast(:,k) .* (g .^ gamma_fast(k)));
end
[Qf, ~] = qr(P_tar_fast, 0);
[Ual, ~, Val] = svd(Qf.' * P_ref_fast, 'econ');
P_tar_fast = Qf * (Ual * Val.');

% SLOW transfer
P_tar_slow = zeros(Ng, r_slow);
for k = 1:r_slow
    P_tar_slow(:,k) = Mmix * (P_ref_slow(:,k) .* (g .^ gamma_slow(k)));
end
[Qs, ~] = qr(P_tar_slow, 0);
[Ual, ~, Val] = svd(Qs.' * P_ref_slow, 'econ');
P_tar_slow = Qs * (Ual * Val.');

% reconstruct
dV_fast = A_ref_fast * P_tar_fast.';    % (T x Ng)
dV_slow = A_ref_slow * P_tar_slow.';    % (T x Ng)

% dip scaling only on FAST part
kappa_det = 0.30;
g_dip_det = (g).^kappa_det;
dV_fast = dV_fast .* repmat(g_dip_det.', T, 1);

% blended deterministic
dV_hat_tar_det = (1-BlendT).*dV_fast + BlendT.*dV_slow;

%% Steady-state operating-point shift 

t_ss_start = 15.0;
beta_ss    = 0.4;
lambda_ss  = 2.0;
tau_ss_sec = 1.0;

idx_ss = find(time >= t_ss_start, 1, 'first'):T;
if isempty(idx_ss)
    idx_ss = round(0.75*T):T;
end

vss_ref = mean(dV_ref(idx_ss,:), 1).';     % (Ng x 1)
v_prior = (g.^beta_ss) .* vss_ref;         % (Ng x 1)

Aeq = eye(Ng) + lambda_ss * Lg;
vss_tar = Aeq \ v_prior;                   % (Ng x 1)

dt = mean(diff(time));
t0 = N0 + 1;
tau_ss = max(round(tau_ss_sec / max(dt,1e-6)), 1);

gate_ss = zeros(T,1);
for t = t0:T
    gate_ss(t) = 1 - exp(-(t - t0)/tau_ss);
end
gate_ss(1:N0) = 0;

for t = 1:T
    if mean(abs(dV_hat_tar_det(t,:)))<1e-4
        gate_ss(t) = 0;
    end
end

% apply shift with reduced gain
ss_gain = 0.5;  % tune: 0~1
dV_hat_tar_det = dV_hat_tar_det + ss_gain * gate_ss*(vss_tar.');

figure;
plot(time, dV_tar_true, 'k'); hold on;
plot(time, dV_hat_tar_det, 'r--')
xlabel('Time (sec)'); ylabel('Voltage (p.u.)'); fontsize(24,"points");

%% Load Bus Voltage Estimation
genvolts = dV_hat_tar_det + V0_tar;
dgenvolts = dV_hat_tar_det;

[ori_prop, H_G, u_g, u_l, s, v] = lv_estm(H, dgenvolts);

lv_prop = ori_prop + Volt3(1,1:29);
lv_real = Volt3(:, 1:29);

%% Voltage Estimation Figures

Err_g_prop = abs(dV_hat_tar_det-dV_tar_true);
mae_g_prop = mean(mean(Err_g_prop));

cmin = 0;
cmax = 0.095;

figure;
imagesc(time, 1:size(Err_g_prop, 2), Err_g_prop')
set(gca,'YDir','normal')
caxis([cmin cmax])
xlabel('Time (sec)'); ylabel('Generator bus #'); fontsize(24,"points");
colormap(jet)
colorbar

Err_prop = abs(lv_real-lv_prop);
mae_prop = mean(mean(Err_prop));

figure;
imagesc(time, 1:size(Err_prop, 2), Err_prop')
set(gca,'YDir','normal')
caxis([cmin cmax])
xlabel('Time (s)'); ylabel('Load bus #'); fontsize(24,"points");
colormap(jet)
colorbar

%% Parameter Settings
cont_size = 500; % Contingency size (Base case: 500 MW)

% Governor and turbine model parameters
Rg = 0.05;            % Governor droop
Tg = 0.5;             % Governor time constant
Trh = 7;              % Reheat time constant
Fh = 0.3;             % Fraction of reheat power

% GFL Coefficients
RGFL = 0.045;
Tconv = 0.02; % Combined with measurement & activation delays
KGFL = 10;

% GFM Coefficients
Ome_b = 2*pi*60;
K_GFM = 3;
Droop_GFM = 0.045;
H_gfm = 4;
D_gfm = 2;

% FFR Coefficients
Rffr = 0.05;
Tffr_dev = 0.02;
Tffr_roc = 0.02;

%% Dynamic Equations: System, SG, and Contingency
% System Dynamics
Num = tf(1, [2*Hsys, Dsys]);
[Anum, Bnum, Cnum, Dnum] = ssdata(Num);

% Generator
Gen = tf(1, [Tg, 1]) * tf([Fh * Trh, 1], [Trh, 1]) * (1 / Rg);
% AGC = tf(5, [1, 0]); 
% Gen1 = parallel(Gen, AGC); % Gen1 can be used instead for SG with AGC
[Agen, Bgen, Cgen, Dgen] = ssdata(Gen);

% Load Contingency
P_m = @(t) 0 * (t < 1) + (-cont_size/Sbase) * (t >= 1); % P_a is 0 until t = 1, then -0.1

DER_size = 600;

btm_load = rand(size(ini_load_v));
btm_load = btm_load/sum(btm_load)*DER_size;
load_p = ini_load_p + btm_load;

v_corr_only_load_bus = lv_prop(:, load_bus);
v_corr_gen_load_bus = genvolts(:, gen_load_bus);
v_load_bus = [v_corr_only_load_bus, v_corr_gen_load_bus];

% Load Contingency
P_m = @(t) 0 * (t < 1) + (-cont_size/Sbase) * (t >= 1); % P_a is 0 until t = 1, then -0.1

% System Dynamics
Num = tf(1, [2*Hsys, Dsys]);
[Anum, Bnum, Cnum, Dnum] = ssdata(Num);

% Generator
Gen = tf(1, [Tg, 1]) * tf([Fh * Trh, 1], [Trh, 1]) * (1 / Rg);
% AGC = tf(5, [1, 0]);
% 
% Gen1 = parallel(Gen, AGC);
[Agen, Bgen, Cgen, Dgen] = ssdata(Gen);

% gfl
GFLsize = GFLbase + DER_size;
GFL = GFLsize/Sbase/RGFL * tf(1, [Tconv, 1]) * tf(1, [Tconv, 1])* tf(KGFL, [1, 0]) / (1+tf(KGFL, [1, 0]));
[AGFL, BGFL, CGFL, DGFL] = ssdata(GFL);

%% Deterministic SFR
load_a = 0.4;
load_b = 0.4;
load_c = 1-load_a-load_b;
loads = Sbase;
load_riv = @(k, riv) sum(load_p.*(load_a.*(riv(k,:)./ini_load_v).^2 + load_b.*(riv(k,:)./ini_load_v) + load_c - 1))/Sbase;
t = 0:dt:T+dt;

x_gen = zeros(size(Agen, 1), N);
x_sim = zeros(size(Anum, 1), N);
x_gfl = zeros(size(AGFL, 1), N);

y_gen = zeros(1,N);
y_gfl = zeros(1,N);
y_sim = zeros(1,N);

lc = zeros(1,N);

for i = 2:N   
    y_all = y_gen(:, i-1)+y_gfl(:, i-1)+lc(:, i-1);
    u = P_m(t(i-1))-y_all;

    [x_sim(:, i), y_sim(:, i)] = runxy(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt);

    % volts
    lc(:,i) = load_riv(i, v_load_bus);
    [x_gen(:, i), y_gen(:, i)] = runxy(x_gen(:, i-1), y_sim(:,i), Agen, Bgen, Cgen, Dgen, dt);
    [x_gfl(:, i), y_gfl(:, i)] = runxy(x_gfl(:, i-1), y_sim(:,i), AGFL, BGFL, CGFL, DGFL, dt);
    y_gfl(:, i) = max(-0.4*GFLsize/Sbase, y_gfl(:,i));
end

freq_det = (1+y_sim)*60;


%% Probabilistic SFR
Nrsim = 10;

% Vg uncertainty parameters
sigma_A0 = 0.6 * (1 + 0.10*d_loc);  % tune: 0.02~0.08
tau_sec  = 0.30;                    % time correlation (sec)
Ltime    = max(5, round(tau_sec / max(dt,1e-6)));

Afast_scale = max(abs(A_ref_fast), [], 1);  Afast_scale = max(Afast_scale, 1e-4);
Aslow_scale = max(abs(A_ref_slow), [], 1);  Aslow_scale = max(Aslow_scale, 1e-4);

sigma_P0   = 0.2 * (1 + 0.05*d_loc);  % tune: 0.05~0.15
beta_diff  = 0.30;
Md = (eye(Ng) + beta_diff*Lg) \ eye(Ng);

Pfast_scale = max(abs(P_tar_fast), [], 1);  Pfast_scale = max(Pfast_scale, 1e-3);
Pslow_scale = max(abs(P_tar_slow), [], 1);  Pslow_scale = max(Pslow_scale, 1e-3);

amp_clip = 3.0;

amp_det = sqrt(mean(dV_hat_tar_det.^2, 2));
amp_det = movmean(amp_det, max(5,round(0.30/dt)));
amp_det_n = amp_det / (max(amp_det)+1e-12);
env_t = min(max(amp_det_n, 0.02), 1.0);
env_t(1:N0) = 0;

slow_factor = 0.4;

rec_prop = zeros([size(data,1), size(lv_real, 2), Nrsim]);
rec_g_prop = zeros([size(data,1), size(V_tar, 2), Nrsim]);

unc_der_s = 1;
unc_lmp_d = 1;
unc_v_d = 1;
unc_la_s = 1;
unc_v = 1;


freq_rec = zeros(size(data,1)-1, Nrsim);
y_gen_rec = zeros(size(data,1)-1, Nrsim);
y_load_rec = zeros(size(data,1)-1, Nrsim);
y_gfl_rec = zeros(size(data,1)-1, Nrsim);

for sim=1:Nrsim

    % BTM & Pure Load
    if unc_der_s == 1
        DER_size = 600*(1+2*(rand()-0.5));        
    else
        DER_size = 600;
    end 

    btm_load = rand(size(ini_load_v));
    btm_load = btm_load/sum(btm_load)*DER_size; % Dirichlet distributed BTM DERs

    load_p = ini_load_p + btm_load;

    if unc_la_s == 1
        load_p = load_p*(1+(rand()-0.5)*0.2);
    end

    % GFL & GFM Inverters
    GFLsize = GFLbase + DER_size; % Large GFL + DER (DERs are also assumed to do GFL ctrl)
    GFMsize = 0; % Size of droop-based GFM
    VSMsize = 0; % Size of VSM GFM

    GFL = GFLsize/Sbase/RGFL * tf(1, [Tconv, 1]) * tf(1, [Tconv, 1])* tf(KGFL, [1, 0]) / (1+tf(KGFL, [1, 0]));
    GFM_droop = GFMsize/Sbase*K_GFM*tf(1, [Tconv, 1])*tf(Ome_b, [1, 0])*tf(1, [Tconv, 1])/(1+K_GFM*Droop_GFM*tf(Ome_b, [1, 0]));
    J_gfm = 2*H_gfm*VSMsize/(Ome_b^2);
    GFM_vsm = VSMsize/Sbase*K_GFM*tf(1, [Tconv, 1])*tf(Ome_b, [1, 0])*tf(1, [Tconv, 1])/(1+tf(1, [J_gfm, D_gfm])*Droop_GFM*tf(Ome_b, [1, 0]));

    [AGFM_droop, BGFM_droop, CGFM_droop, DGFM_droop] = ssdata(GFM_droop);
    [AGFM_vsm, BGFM_vsm, CGFM_vsm, DGFM_vsm] = ssdata(GFM_vsm);
    [AGFL, BGFL, CGFL, DGFL] = ssdata(GFL);

    % FFR Resources
    FFR_dev_size = 0;
    FFR_roc_size = 0;
    FFR_dev = tf(1, [Tffr_dev, 1]) * tf(1, [Tconv, 1]) * (1 / Rffr) * (FFR_dev_size)/Sbase;
    [Affrdev, Bffrdev, Cffrdev, Dffrdev] = ssdata(FFR_dev);

    FFR_roc = tf([1,0], 1) * tf(1, [Tffr_dev, 1]) * tf(1, [Tconv, 1]) * (1 / Rffr) * FFR_roc_size/Sbase;
    [Affrroc, Bffrroc, Cffrroc, Dffrroc] = ssdata(FFR_roc);

    % Voltage probabilistic generation
    if unc_v == 1
        Et_fast = randn(T, r_fast);
        Et_fast = movmean(Et_fast, Ltime, 1);
        Et_fast(1:N0,:) = 0;
    
        Et_slow = randn(T, r_slow);
        Et_slow = movmean(Et_slow, Ltime, 1);
        Et_slow(1:N0,:) = 0;
    
        % apply envelope (vanish in steady-state)
        Et_fast = Et_fast .* (sigma_A0 * (env_t * Afast_scale));
        Et_slow = Et_slow .* (slow_factor*sigma_A0 * (env_t * Aslow_scale));
    
        % clip
        Et_fast = min(max(Et_fast, -amp_clip*sigma_A0*repmat(Afast_scale,T,1)), ...
                              amp_clip*sigma_A0*repmat(Afast_scale,T,1));
        Et_slow = min(max(Et_slow, -amp_clip*slow_factor*sigma_A0*repmat(Aslow_scale,T,1)), ...
                              amp_clip*slow_factor*sigma_A0*repmat(Aslow_scale,T,1));
    
        A_fast_sto = A_ref_fast + Et_fast;
        A_slow_sto = A_ref_slow + Et_slow;
    
        %% (B) Spatial noise on P_fast and P_slow (graph-smoothed)
        Ep_fast = randn(Ng, r_fast);
        Ep_fast = Md * Ep_fast;
        Ep_fast = Ep_fast .* (sigma_P0 * repmat(Pfast_scale, Ng, 1));
        Ep_fast = min(max(Ep_fast, -amp_clip*sigma_P0*repmat(Pfast_scale,Ng,1)), ...
                              amp_clip*sigma_P0*repmat(Pfast_scale,Ng,1));
        P_fast_sto = P_tar_fast + Ep_fast;
    
        Ep_slow = randn(Ng, r_slow);
        Ep_slow = Md * Ep_slow;
        Ep_slow = Ep_slow .* (sigma_P0 * repmat(Pslow_scale, Ng, 1));
        Ep_slow = min(max(Ep_slow, -amp_clip*sigma_P0*repmat(Pslow_scale,Ng,1)), ...
                              amp_clip*sigma_P0*repmat(Pslow_scale,Ng,1));
        P_slow_sto = P_tar_slow + Ep_slow;
    
        % Align each noisy P to its deterministic subspace (keeps modes stable)
        [Qf2,~] = qr(P_fast_sto,0);
        [Ual2,~,Val2] = svd(Qf2.'*P_tar_fast,'econ');
        P_fast_sto = Qf2*(Ual2*Val2.');
    
        [Qs2,~] = qr(P_slow_sto,0);
        [Ual2,~,Val2] = svd(Qs2.'*P_tar_slow,'econ');
        P_slow_sto = Qs2*(Ual2*Val2.');
    
        %% (C) Reconstruct fast/slow and blend
        dV_fast_sto = A_fast_sto * P_fast_sto.';
        dV_slow_sto = A_slow_sto * P_slow_sto.';
    
        % dip scaling only on fast (same as deterministic)
        dV_fast_sto = dV_fast_sto .* repmat(g_dip_det.', T, 1);
    
        dV_hat_tar = (1-BlendT).*dV_fast_sto + BlendT.*dV_slow_sto;
    
        % apply SS shift (reduced)
        dV_hat_tar = dV_hat_tar + ss_gain * gate_ss*(vss_tar.');
    
        % strict no randomness pre-disturbance: match deterministic baseline
        dV_hat_tar(1:N0,:) = dV_hat_tar_det(1:N0,:);
        
        pro_prop = propvoltage(dV_hat_tar_det, H_G, u_g, u_l, s, v) + Volt3(1,1:29);
    else
        pro_prop = lv_prop;
    end

    pro_lbus = pro_prop(:, load_bus);
    v_load_bus_prob = [pro_lbus, v_corr_gen_load_bus];

    % Voltage Dependency of Load (ZIP Load Model)
    if unc_lmp_d == 1
        load_a = 0.4+0.2*(rand()-0.5); % Coefficient Z
        load_b = 0.4+0.2*(rand()-0.5); % Coefficient I
        load_c = 1-load_a-load_b;      % Coefficient P
    else
        load_a = 0.4; % Coefficient Z
        load_b = 0.4; % Coefficient I
        load_c = 1-load_a-load_b;      % Coefficient P
    end
    load_riv = @(k, riv) sum(load_p.*(load_a.*(riv(k,:)./ini_load_v).^2 + load_b.*(riv(k,:)./ini_load_v) + load_c - 1))/(sum(load_p)+cont_size);

    t = 0:dt:T+dt;

    x_gen = zeros(size(Agen, 1), N);
    x_sim = zeros(size(Anum, 1), N);
    x_gfl = zeros(size(AGFL, 1), N);
    x_gfm = zeros(size(AGFM_droop, 1), N);
    x_vsm = zeros(size(AGFM_vsm, 1), N);
    x_ffrdev = zeros(size(Affrdev, 1), N);
    x_ffrroc = zeros(size(Affrroc, 1), N);
    trip_size = zeros(1,N);

    y_gen = zeros(1,N);
    y_gfl = zeros(1,N);
    y_gfm = zeros(1,N);
    y_vsm = zeros(1,N);
    y_sim = zeros(1,N);
    y_ffrdev = zeros(1,N);
    y_ffrroc = zeros(1,N);

    lc = zeros(1,N);
    lc2 = zeros(1,N);
    for i = 2:N   
        % 
        y_all = y_gen(:, i-1)+y_gfl(:, i-1)+lc(:, i-1)+y_gfm(:, i-1)+y_vsm(:, i-1)+y_ffrdev(:, i-1)+y_ffrroc(:, i-1);
        u = P_m(t(i-1))-y_all;

        % To include unknown uncertainty, use runsde; Else, use runxy
        % [x_sim(:, i), y_sim(:, i)] = runxy(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt);
        [x_sim(:, i), y_sim(:, i)] = runsde(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt, t(i));

        % Voltages & DER Trips
        if unc_v_d == 1
            v_load_bus_prob(i,:) = v_load_bus_prob(i,:)-0.1*v_load_bus_prob(i,:).*rand(1,19).*(trip_size(i));
        end

        lc(:,i) = load_riv(i, v_load_bus_prob);
        trip_signal(i,:) = DER_vrt(size(v_load_bus_prob,2), t(i), dt, v_load_bus_prob(i,:)); % All DERs are assumed to be legacy DER in here
        trip_size(i) = sum(btm_load.*(trip_signal(i,:)>0));

        [x_gen(:, i), y_gen(:, i)] = runxy(x_gen(:, i-1), y_sim(:,i), Agen, Bgen, Cgen, Dgen, dt);
        [x_gfl(:, i), y_gfl(:, i)] = runxy(x_gfl(:, i-1), y_sim(:,i), AGFL, BGFL, CGFL, DGFL, dt);
        [x_gfm(:, i), y_gfm(:, i)] = runxy(x_gfm(:, i-1), y_sim(:,i), AGFM_droop, BGFM_droop, CGFM_droop, DGFM_droop, dt);
        [x_vsm(:, i), y_vsm(:, i)] = runxy(x_vsm(:, i-1), y_sim(:,i), AGFM_vsm, BGFM_vsm, CGFM_vsm, DGFM_vsm, dt);
        [x_ffrdev(:,i), y_ffrdev(:,i)] = runxy(x_ffrdev(:, i-1), y_sim(:,i), Affrdev, Bffrdev, Cffrdev, Dffrdev, dt);
        [x_ffrroc(:,i), y_ffrroc(:,i)] = runxy(x_ffrroc(:, i-1), y_sim(:,i), Affrroc, Bffrroc, Cffrroc, Dffrroc, dt);

        y_gfl(:, i) = max(-0.4*GFLsize/Sbase, y_gfl(:,i));
        y_gfl(:, i) = y_gfl(:,i)*(GFLsize-trip_size(i))/GFLsize;
        y_gfm(:, i) = max(-0.4*GFMsize/Sbase, y_gfm(:,i));
        y_vsm(:, i) = max(-0.4*VSMsize/Sbase, y_vsm(:,i));

    end
    freq = (1+y_sim)*60;
    min_freq(sim) = min(freq);

    % Recording
    freq_rec(:, sim) = freq;
    y_gen_rec(:,sim) = -y_gen;
    y_gfl_rec(:,sim) = -y_gfl;
    y_load_rec(:,sim) = -lc;
    rec_prop(:,:,sim) = pro_prop;
    rec_g_prop(:,:,sim) = dV_hat_tar + V0_tar;
end

[mean_freq, ci_z, ci_m] = find_ci(freq_rec);
[mean_gen, ci_g, ci_g2] = find_ci(y_gen_rec);
[mean_gfl, ci_i, ci_i2] = find_ci(y_gfl_rec);
[mean_load, ci_l, ci_l2] = find_ci(y_load_rec);
[mean_v, ci_v, ci_v2] = find_ci(rec_prop);
[mean_vg, ci_vg, ci_vg2] = find_ci(rec_g_prop);

convt = time(1:2402);

lind = [3, 23];
gind = [1, 10];

figure;
for i=1:size(lind,2)
    lvi = lind(i);
    specvs = reshape(rec_prop(:,lvi,:), [size(rec_prop,1), Nrsim]);
    [mean_Vs, ci_zs, ci_ms] = find_ci_robust(specvs);
    fill([time; flipud(time)], [ci_ms(:,1); flipud(ci_ms(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'r', 'EdgeAlpha', 0.15); hold on;
    fill([time; flipud(time)], [ci_zs(:,1); flipud(ci_zs(:,2))], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.3);
end
plot(time, Volt3(:,lind), 'Color', 'k', 'LineWidth', 2); hold on;
plot(time, lv_prop(:,lind), 'LineStyle', ':', 'LineWidth', 2, 'Color', 'r'); title('Load bus');
xlabel('Time (s)'); ylabel('Voltage (p.u.)'); fontsize(24,"points");

figure;
for i=1:size(gind,2)
    lvi = gind(i);
    specvs = reshape(rec_g_prop(:,lvi,:), [size(rec_g_prop,1), Nrsim]);
    [mean_Vs, ci_zs, ci_ms] = find_ci_robust(specvs);
    fill([time; flipud(time)], [ci_ms(:,1); flipud(ci_ms(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'r', 'EdgeAlpha', 0.15); hold on;
    fill([time; flipud(time)], [ci_zs(:,1); flipud(ci_zs(:,2))], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.3);
end
plot(time, V_tar(:,gind), 'Color', 'k', 'LineWidth', 2); hold on;
plot(time, genvolts(:,gind), 'LineStyle', ':', 'LineWidth', 2, 'Color', 'r');  title('Generator bus');
xlabel('Time (s)'); ylabel('Voltage (p.u.)'); fontsize(24,"points");

figure;
fill([convt; flipud(convt)], [ci_m(:,1); flipud(ci_m(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'r', 'EdgeAlpha', 0.15); hold on;
fill([convt; flipud(convt)], [ci_z(:,1); flipud(ci_z(:,2))], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.3);
plot(convt, (freqsys(2:end)+1)*60, 'Color', 'k', 'LineWidth', 2); hold on;
plot(convt, freq_det, 'LineStyle', ':', 'LineWidth', 2, 'Color', 'r');
xlabel('Time (s)'); ylabel('Frequency (Hz)'); fontsize(24,"points");

figure;
fill([convt; flipud(convt)], [ci_g2(:,1); flipud(ci_g2(:,2))], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'EdgeAlpha', 0.15); hold on;
fill([convt; flipud(convt)], [ci_g(:,1); flipud(ci_g(:,2))], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'EdgeAlpha', 0.3);
fill([convt; flipud(convt)], [ci_i2(:,1); flipud(ci_i2(:,2))], 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'b', 'EdgeAlpha', 0.15); hold on;
fill([convt; flipud(convt)], [ci_i(:,1); flipud(ci_i(:,2))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'EdgeAlpha', 0.3);
fill([convt; flipud(convt)], [ci_l2(:,1); flipud(ci_l2(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'r', 'EdgeAlpha', 0.15); hold on;
fill([convt; flipud(convt)], [ci_l(:,1); flipud(ci_l(:,2))], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.3);
plot(convt, mean_gen, 'k', 'LineWidth', 3);
plot(convt, mean_gfl, 'b', 'LineWidth', 3, 'LineStyle', '--');
plot(convt, mean_load, 'r', 'LineWidth', 3, 'LineStyle', '-.');
% legend('Generator', 'GFL', 'Load')
xlabel('Time (s)'); ylabel('Normalized P'); fontsize(24,"points")

[mu,s,muci,sci] = normfit(min_freq);
% gamma = fitdist(min_freq', 'Gamma');

x = min(min_freq):0.001:max(min_freq);
y = normpdf(x,mu,s);
% y_gamma = gampdf(x,gamma.a,gamma.b);

pd = fitdist(min_freq', 'Normal');
ci = paramci(pd,'Alpha',.01);   

figure;
h = histogram(min_freq, 30, 'Orientation', 'horizontal','Normalization','pdf', 'facecolor', 'r', 'FaceAlpha', 0.5);
hold on
plot(y,x, 'r','LineWidth', 2)
yline(min((freqsys(2:end)+1)*60), 'k','LineWidth', 2)
ylabel('Frequency Nadir'); xlabel('# of bins'); fontsize(24,"points");

figure;
h = histogram(min_freq, 30, 'Orientation', 'vertical','Normalization','pdf', 'facecolor', 'r', 'FaceAlpha', 0.5);
hold on
plot(x,y, 'LineWidth', 2)
xline(min((freqsys(2:end)+1)*60), 'k','LineWidth', 2)
ylabel('Frequency Nadir'); xlabel('# of bins'); fontsize(24,"points");

%% Local Functions

function [x2, y2] = runxy(x, u, A, B, C, D, dt)
   drift = A*x + B*u;
   x2 = x + drift*dt;
   y2 = C*x2 + D*u;
end

function [x2, y2] = runsde(x, u, A, B, C, D, dt,t)
   drift = A*x + B*u;
   diffusion_f = (0.3*x*(rand()-0.5))*(t>1);
   noise_f = sqrt(dt)*rand(size(drift));
   diffusion_p = (0.15*u*(rand()-0.5))*(t>1);
   noise_p = sqrt(dt)*rand(size(drift));
   x2 = x + drift*dt+ diffusion_p*noise_p;
   y2 = C*x2 + D*u + diffusion_f*noise_f;
end

function proposed = propvoltage(dgenvolts, H_G, u_g, u_l, s, v)
    sigRow = 0.03;        % row scaling
    sigCol = 0.03;        % col scaling
    for i=1:size(dgenvolts,1)
        dV = dgenvolts(i,:);
        dQest(i,:) = pinv(H_G)*dV';
        hat_VG = H_G*dQest(i,:)';
        residual(i,:) = dV-hat_VG';
    
        x = v'*dQest(i,:)';
        y = u_g'*dV';
    
        e_t = y - s*x;
        ds_t = e_t./x;
        s2_base = s+diag(ds_t);
        H_Gt = u_g*s2_base*v';
        H_Lt = u_l*s2_base*v';
        S_LGt0 = H_Lt*pinv(H_Gt);
        
        nL = size(S_LGt0,1);
        nG = size(S_LGt0,2);

        dL = exp(sigRow*randn(nL,1) - 0.5*sigRow^2);  % mean=1
        dG = exp(sigCol*randn(nG,1) - 0.5*sigCol^2);  % mean=1

        S_LGt = (dL .* S_LGt0) .* (dG.');  % same as diag(dL)*S0*diag(dG)
        proposed(i,:) = S_LGt * dV';
    end
end

function vrt_out = DER_vrt(N, t, dt, loadvolt)
    persistent v_buffer
    persistent der_trip_rec

    buffer_size_1 = round(0.5/dt); % 0.5 sec under 0.88 p.u.

    if isempty(v_buffer)
        v_buffer = zeros(N,1);
    end

    if isempty(der_trip_rec)
        der_trip_rec = zeros(N,1);
    end
    if t<0.1
        v_buffer = zeros(N,1);
        der_trip_rec = zeros(N,1);
    end

    for b=1:size(loadvolt, 2)
        if der_trip_rec(b)==1
            der_trip_rec(b)=1;
        else
            uv1 = (loadvolt(b)<0.88);

            if uv1
                v_buffer(b,1) = v_buffer(b,1)+1;
            else
                v_buffer(b,:) = 0;
            end

            if (v_buffer(b,1)>buffer_size_1)
                der_trip_rec(b) = 1;
            end
        end
    end
    vrt_out = der_trip_rec;
end

function [dist_matrix, net] = build_dist_matrix_from_psse_raw(rawFile, weightMode)
if nargin < 2, weightMode = 'X'; end

lines = readlines(rawFile);
lines = string(lines);

iBusBeg = find(contains(lines, "BEGIN BUS DATA"), 1, 'first');
iBusEnd = find(contains(lines, "END OF BUS DATA"), 1, 'first');
iBrBeg  = find(contains(lines, "BEGIN BRANCH DATA"), 1, 'first');
iBrEnd  = find(contains(lines, "END OF BRANCH DATA"), 1, 'first');
iTxBeg  = find(contains(lines, "BEGIN TRANSFORMER DATA"), 1, 'first');
iTxEnd  = find(contains(lines, "END OF TRANSFORMER DATA"), 1, 'first');

if isempty(iBusBeg) || isempty(iBusEnd) || isempty(iBrBeg) || isempty(iBrEnd)
    error("Cannot locate BUS/BRANCH sections in RAW file.");
end

busIds = parse_bus_ids(lines(iBusBeg:iBusEnd));
Nb = max(busIds);
if Nb <= 0, error("Failed to infer Nb."); end

[br_from, br_to, br_r, br_x] = parse_branch_section(lines(iBrBeg:iBrEnd));

tx_from = []; tx_to = []; tx_r = []; tx_x = [];
if ~isempty(iTxBeg) && ~isempty(iTxEnd)
    [tx_from, tx_to, tx_r, tx_x] = parse_transformer_section(lines(iTxBeg:iTxEnd));
end

from = [br_from; tx_from];
to   = [br_to;   tx_to];
r    = [br_r;    tx_r];
x    = [br_x;    tx_x];

ok = from>=1 & from<=Nb & to>=1 & to<=Nb & isfinite(x) & x>0;
from = from(ok); to = to(ok); r = r(ok); x = x(ok);

net.Nb = Nb; net.from = from; net.to = to; net.r = r; net.x = x;

switch upper(weightMode)
    case 'X'
        w = x;
    case 'ZMAG'
        w = sqrt(r.^2 + x.^2);
    otherwise
        error("weightMode must be 'X' or 'ZMAG'.");
end

G = graph(from, to, w, Nb);
dist_matrix = distances(G);
dist_matrix(~isfinite(dist_matrix)) = 1e6;
end

function busIds = parse_bus_ids(lines)
busIds = [];
for k = 1:numel(lines)
    s = strtrim(lines(k));
    if s == "" || startsWith(s, "@") || startsWith(s, "0 /")
        continue;
    end
    toks = split(s, ",");
    if numel(toks) < 2, continue; end
    id = str2double(strtrim(toks(1)));
    if isfinite(id) && id > 0
        busIds(end+1,1) = id;
    end
end
end

function [from, to, r, x] = parse_branch_section(lines)
from = []; to = []; r = []; x = [];
for k = 1:numel(lines)
    s = strtrim(lines(k));
    if s == "" || startsWith(s, "@") || startsWith(s, "0 /")
        continue;
    end
    toks = split(s, ",");
    if numel(toks) < 5, continue; end
    I = str2double(strtrim(toks(1)));
    J = str2double(strtrim(toks(2)));
    R = str2double(strtrim(toks(4)));
    X = str2double(strtrim(toks(5)));
    if isfinite(I) && isfinite(J) && isfinite(X) && X > 0
        from(end+1,1) = I;
        to(end+1,1)   = J;
        r(end+1,1)    = max(R, 0);
        x(end+1,1)    = X;
    end
end
end

function [from, to, r, x] = parse_transformer_section(lines)
from = []; to = []; r = []; x = [];
k = 1;
while k <= numel(lines)
    s1 = strtrim(lines(k));
    if s1 == "" || startsWith(s1, "@") || startsWith(s1, "0 /")
        k = k + 1;
        continue;
    end
    toks1 = split(s1, ",");
    if numel(toks1) < 4
        k = k + 1;
        continue;
    end
    I = str2double(strtrim(toks1(1)));
    J = str2double(strtrim(toks1(2)));

    if k+1 > numel(lines), break; end
    s2 = strtrim(lines(k+1));
    toks2 = split(s2, ",");
    if numel(toks2) < 2
        k = k + 1;
        continue;
    end
    R12 = str2double(strtrim(toks2(1)));
    X12 = str2double(strtrim(toks2(2)));

    if isfinite(I) && isfinite(J) && isfinite(X12) && X12 > 0
        from(end+1,1) = I;
        to(end+1,1)   = J;
        r(end+1,1)    = max(R12, 0);
        x(end+1,1)    = X12;
    end

    k = k + 3;
end
end

function [mean_x, ci_z, ci_m] = find_ci_robust(x_rec)
% x_rec: (T x Nsamples)
mean_x = mean(x_rec, 2, 'omitnan');
ci_z = zeros(size(x_rec,1), 2);
ci_m = zeros(size(x_rec,1), 2);

for tt=1:size(x_rec,1)
    A = x_rec(tt, :);
    A = A(isfinite(A));
    if isempty(A)
        ci_z(tt,:) = [NaN NaN];
        ci_m(tt,:) = [NaN NaN];
    else
        ci_z(tt,:) = [prctile(A, 5), prctile(A, 95)];
        ci_m(tt,:) = [min(A), max(A)];
    end
end
end

function [mean_freq, ci_z, ci_m] = find_ci(freq_rec)
    mean_freq = mean(freq_rec, 2);
    ci_z = zeros(size(freq_rec,1), 2);
    ci_m = zeros(size(freq_rec,1), 2);
    
    for tt=1:size(freq_rec,1)
        A = freq_rec(tt, :);
        ci_z(tt,:) = [prctile(A, 5), prctile(A, 95)];
        ci_m(tt,:) = [min(A), max(A)];
    end
end
function [ori_prop, H_G, u_g, u_l, s, v] = lv_estm(H, dgenvolts)
    % Proposed Method
    H_L = H(1:29, :);
    H_G = H(30:end, :);
    
    [u, s, v] = svd(H);
    M = 10;
    u_l = u(1:29, 1:M);
    u_g = u(30:end, 1:M);
    s = s(1:M, :);
    
    for i=1:size(dgenvolts,1)
        dV = dgenvolts(i,:);
        dQest(i,:) = pinv(H_G)*dV';
        hat_VG = H_G*dQest(i,:)';
        residual(i,:) = dV-hat_VG';
    
        x = v'*dQest(i,:)';
        y = u_g'*dV';
    
        e_t = y - s*x;
        ds_t = e_t./x;
        s2_base = s+diag(ds_t);
        H_Gt = u_g*s2_base*v';
        H_Lt = u_l*s2_base*v';
        S_LGt0 = H_Lt*pinv(H_Gt);
        ori_prop(i,:) = S_LGt0*dV';
    end
    
    ori_prop(1,:)=0;
end