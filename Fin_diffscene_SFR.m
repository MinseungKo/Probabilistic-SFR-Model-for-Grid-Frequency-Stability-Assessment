clear all
close all
clc

%% Import Data
tic;
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
% Conventional load-bus voltage recovery for the target scenario.
lv_conv = convV_from_dgen(dgenvolts, vsens) + Volt3(1,1:29);

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
v_load_bus_conv = [lv_conv(:, load_bus), v_corr_gen_load_bus];

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


%% Probabilistic SFR and Uncertainty Sensitivity
Nrsim = 500;
unc_scale_list = [0, 0.5, 1.0, 1.5, 2.0];
base_unc_scale = 1.0;
sens_prob_diff = table();
base_prob_store = struct();


for iu = 1:numel(unc_scale_list)
    tic;
    unc_scale = unc_scale_list(iu);
    fprintf('Running different-scenario probabilistic SFR with uncertainty scale = %.2f\n', unc_scale);

    % Vg uncertainty parameters
    sigma_A0 = unc_scale * 0.6 * (1 + 0.10*d_loc);  % tune: 0.02~0.08
    tau_sec  = 0.30;                    % time correlation (sec)
    Ltime    = max(5, round(tau_sec / max(dt,1e-6)));

    Afast_scale = max(abs(A_ref_fast), [], 1);  Afast_scale = max(Afast_scale, 1e-4);
    Aslow_scale = max(abs(A_ref_slow), [], 1);  Aslow_scale = max(Aslow_scale, 1e-4);

    sigma_P0   = unc_scale * 0.2 * (1 + 0.05*d_loc);  % tune: 0.05~0.15
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
            DER_size = 600*(1 + unc_scale*2*(rand()-0.5));        
        else
            DER_size = 600;
        end 

        if unc_scale == 0
            btm_load = ini_load_p/sum(ini_load_p)*DER_size;
        else
            btm_load = rand(size(ini_load_v));
            btm_load = btm_load/sum(btm_load)*DER_size; % Dirichlet distributed BTM DERs
        end

        load_p = ini_load_p + btm_load;

        if unc_la_s == 1
            load_p = load_p*(1 + unc_scale*(rand()-0.5)*0.2);
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

            pro_prop = propvoltage(dV_hat_tar, H_G, u_g, u_l, s, v, unc_scale) + Volt3(1,1:29);
        else
            pro_prop = lv_prop;
        end

        pro_lbus = pro_prop(:, load_bus);
        v_load_bus_prob = [pro_lbus, v_corr_gen_load_bus];

        % Voltage Dependency of Load (ZIP Load Model)
        if unc_lmp_d == 1
            load_a = 0.4 + unc_scale*0.2*(rand()-0.5); % Coefficient Z
            load_b = 0.4 + unc_scale*0.2*(rand()-0.5); % Coefficient I
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
            [x_sim(:, i), y_sim(:, i)] = runsde(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt, t(i), unc_scale);

            % Voltages & DER Trips
            if unc_v_d == 1
                v_load_bus_prob(i,:) = v_load_bus_prob(i,:) - unc_scale*0.1*v_load_bus_prob(i,:).*rand(1,19).*(trip_size(i));
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
    ref_freq_sens = (freqsys3(2:end) + 1) * 60;
    time_sens = time(:);
    Lsens = min([numel(ref_freq_sens), size(freq_rec,1), numel(time_sens)]);
    Tprob_i = calc_prob_metrics(ref_freq_sens(1:Lsens), freq_rec(1:Lsens,:), time_sens(1:Lsens), "Different-Uncertainty sensitivity", 59.5);
    Tprob_i.UncScale = unc_scale;
    Tprob_i = movevars(Tprob_i, 'UncScale', 'After', 'Case');
    sens_prob_diff = [sens_prob_diff; Tprob_i];

    if abs(unc_scale - base_unc_scale) < 1e-12
        base_prob_store.freq_rec = freq_rec;
        base_prob_store.y_gen_rec = y_gen_rec;
        base_prob_store.y_gfl_rec = y_gfl_rec;
        base_prob_store.y_load_rec = y_load_rec;
        base_prob_store.rec_prop = rec_prop;
        base_prob_store.rec_g_prop = rec_g_prop;
        base_prob_store.min_freq = min_freq;
    end

    toc
end


freq_rec = base_prob_store.freq_rec;
y_gen_rec = base_prob_store.y_gen_rec;
y_gfl_rec = base_prob_store.y_gfl_rec;
y_load_rec = base_prob_store.y_load_rec;
rec_prop = base_prob_store.rec_prop;
rec_g_prop = base_prob_store.rec_g_prop;
min_freq = base_prob_store.min_freq;



%% Probabilistic Baseline Comparison: Different Scenario
[freq_rec_novolt, freq_rec_conv, freq_rec_fullsde] = run_probabilistic_baseline_comparison( ...
    Nrsim, N, t, dt, Sbase, cont_size, P_m, ...
    Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
    GFLbase, RGFL, Tconv, KGFL, ...
    ini_load_v, ini_load_p, 600, v_load_bus_conv, v_load_bus, base_unc_scale);

[mean_freq, ci_z, ci_m] = find_ci(freq_rec);
[mean_gen, ci_g, ci_g2] = find_ci(y_gen_rec);
[mean_gfl, ci_i, ci_i2] = find_ci(y_gfl_rec);
[mean_load, ci_l, ci_l2] = find_ci(y_load_rec);
[mean_v, ci_v, ci_v2] = find_ci(rec_prop);
[mean_vg, ci_vg, ci_vg2] = find_ci(rec_g_prop);

%% Quantitative Metrics: Different Scenario
f_lim = 59.5;

ref_freq_diff = (freqsys3(2:end) + 1) * 60;
time_metric_diff = time(:);

Lmet = min([numel(ref_freq_diff), numel(freq_det), size(freq_rec,1), numel(time_metric_diff)]);
ref_freq_diff = ref_freq_diff(1:Lmet);
time_metric_diff = time_metric_diff(1:Lmet);

freq_det_m = freq_det(1:Lmet);
freq_rec_m = freq_rec(1:Lmet,:);

Tdet_diff = calc_det_metrics(ref_freq_diff, freq_det_m, time_metric_diff, "Different-Proposed deterministic SFR", f_lim);

Vg_ref_diff = V_tar;
Vg_hat_diff = genvolts;
Lg_metric = min([size(Vg_ref_diff,1), size(Vg_hat_diff,1)]);
Vg_ref_diff = Vg_ref_diff(1:Lg_metric,:);
Vg_hat_diff = Vg_hat_diff(1:Lg_metric,:);

gen_voltage_metrics_diff = table( ...
    "Different-Generator-bus voltage recovery", ...
    mean(abs(Vg_hat_diff - Vg_ref_diff), "all"), ...
    sqrt(mean((Vg_hat_diff - Vg_ref_diff).^2, "all")), ...
    max(abs(Vg_hat_diff - Vg_ref_diff), [], "all"), ...
    'VariableNames', {'Case', 'VoltageMAE_pu', 'VoltageRMSE_pu', 'VoltageMaxAE_pu'} ...
);

Vl_ref_diff = lv_real;
Vl_hat_diff = lv_prop;
Ll_metric = min([size(Vl_ref_diff,1), size(Vl_hat_diff,1)]);
Vl_ref_diff = Vl_ref_diff(1:Ll_metric,:);
Vl_hat_diff = Vl_hat_diff(1:Ll_metric,:);

load_voltage_metrics_diff = table( ...
    "Different-Load-bus voltage recovery", ...
    mean(abs(Vl_hat_diff - Vl_ref_diff), "all"), ...
    sqrt(mean((Vl_hat_diff - Vl_ref_diff).^2, "all")), ...
    max(abs(Vl_hat_diff - Vl_ref_diff), [], "all"), ...
    'VariableNames', {'Case', 'VoltageMAE_pu', 'VoltageRMSE_pu', 'VoltageMaxAE_pu'} ...
);

Vl_conv_diff = lv_conv(1:Ll_metric,:);

voltage_det_metrics_diff = [
    calc_voltage_det_metrics(Vg_ref_diff, Vg_hat_diff, time(1:Lg_metric), "Different-Proposed generator-bus voltage recovery");
    calc_voltage_det_metrics(Vl_ref_diff, Vl_hat_diff, time(1:Ll_metric), "Different-Proposed load-bus voltage recovery");
    calc_voltage_det_metrics(Vl_ref_diff, Vl_conv_diff, time(1:Ll_metric), "Different-Conventional load-bus voltage recovery")
];

Tprob_diff = calc_prob_metrics(ref_freq_diff, freq_rec_m, time_metric_diff, "Different-Proposed probabilistic SFR", f_lim);


freq_rec_novolt_m = freq_rec_novolt(1:Lmet,:);
freq_rec_conv_m = freq_rec_conv(1:Lmet,:);
freq_rec_fullsde_m = freq_rec_fullsde(1:Lmet,:);

Tprob_compare_diff = [
    Tprob_diff;
    calc_prob_metrics(ref_freq_diff, freq_rec_novolt_m, time_metric_diff, "Different-Baseline SFR without voltage recovery", f_lim);
    calc_prob_metrics(ref_freq_diff, freq_rec_conv_m, time_metric_diff, "Different-Baseline SFR with conventional voltage recovery", f_lim);
    calc_prob_metrics(ref_freq_diff, freq_rec_fullsde_m, time_metric_diff, "Different-Baseline full-SDE SFR without explicit uncertainty modeling", f_lim)
];

rec_prop_m = rec_prop(1:min(size(rec_prop,1), size(lv_real,1)),:,:);
Vref_load_prob = lv_real(1:size(rec_prop_m,1),:);

qVl025 = prctile(rec_prop_m, 2.5, 3);
qVl975 = prctile(rec_prop_m, 97.5, 3);

Vl_PICP_95 = mean((Vref_load_prob >= qVl025) & (Vref_load_prob <= qVl975), "all");
Vl_PINAW_95 = mean(qVl975 - qVl025, "all") / (max(Vref_load_prob, [], "all") - min(Vref_load_prob, [], "all") + eps);

voltage_prob_metrics_diff_load = table( ...
    "Different-Probabilistic load-bus voltage recovery", Vl_PICP_95, Vl_PINAW_95, ...
    'VariableNames', {'Case', 'VoltagePICP_95', 'VoltagePINAW_95'} ...
);

rec_g_prop_m = rec_g_prop(1:min(size(rec_g_prop,1), size(V_tar,1)),:,:);
Vref_gen_prob = V_tar(1:size(rec_g_prop_m,1),:);

qVg025 = prctile(rec_g_prop_m, 2.5, 3);
qVg975 = prctile(rec_g_prop_m, 97.5, 3);

Vg_PICP_95 = mean((Vref_gen_prob >= qVg025) & (Vref_gen_prob <= qVg975), "all");
Vg_PINAW_95 = mean(qVg975 - qVg025, "all") / (max(Vref_gen_prob, [], "all") - min(Vref_gen_prob, [], "all") + eps);

voltage_prob_metrics_diff_gen = table( ...
    "Different-Probabilistic generator-bus voltage recovery", Vg_PICP_95, Vg_PINAW_95, ...
    'VariableNames', {'Case', 'VoltagePICP_95', 'VoltagePINAW_95'} ...
);

print_metric_table(Tdet_diff, "Deterministic frequency metrics: different scenario");
print_metric_table(gen_voltage_metrics_diff, "Generator-bus voltage recovery metrics: different scenario");
print_metric_table(load_voltage_metrics_diff, "Load-bus voltage recovery metrics: different scenario");
print_metric_table(voltage_det_metrics_diff, "Expanded deterministic voltage recovery metrics: different scenario");
print_metric_table(Tprob_diff, "Probabilistic frequency metrics: different scenario");
print_metric_table(Tprob_compare_diff, "Probabilistic frequency metric comparison: different scenario");
print_metric_table(voltage_prob_metrics_diff_gen, "Probabilistic generator-bus voltage coverage metrics: different scenario");
print_metric_table(voltage_prob_metrics_diff_load, "Probabilistic load-bus voltage coverage metrics: different scenario");
print_metric_table(sens_prob_diff, "Uncertainty sensitivity metrics: different scenario");
% 
safe_writetable(Tdet_diff, "metrics_diff_deterministic_frequency.csv");
safe_writetable(gen_voltage_metrics_diff, "metrics_diff_generator_voltage.csv");
safe_writetable(load_voltage_metrics_diff, "metrics_diff_load_voltage.csv");
safe_writetable(voltage_det_metrics_diff, "metrics_diff_voltage_expanded.csv");
safe_writetable(Tprob_diff, "metrics_diff_probabilistic_frequency.csv");
safe_writetable(Tprob_compare_diff, "metrics_diff_probabilistic_frequency_comparison.csv");
safe_writetable(voltage_prob_metrics_diff_gen, "metrics_diff_probabilistic_generator_voltage.csv");
safe_writetable(voltage_prob_metrics_diff_load, "metrics_diff_probabilistic_load_voltage.csv");
safe_writetable(sens_prob_diff, "metrics_diff_uncertainty_sensitivity.csv");

figure;
plot(sens_prob_diff.UncScale, sens_prob_diff.NadirQ05_Hz, '-o', 'LineWidth', 2); hold on;
plot(sens_prob_diff.UncScale, sens_prob_diff.NadirQ50_Hz, '-s', 'LineWidth', 2);
plot(sens_prob_diff.UncScale, sens_prob_diff.NadirQ95_Hz, '-^', 'LineWidth', 2);
yline(f_lim, 'k--', 'LineWidth', 2);
grid on;
xlabel('Uncertainty scale'); ylabel('Frequency nadir (Hz)');
legend('5% nadir', 'Median nadir', '95% nadir', 'Security limit', 'Location', 'best');
fontsize(20, "points");

figure;
yyaxis left
plot(sens_prob_diff.UncScale, sens_prob_diff.PICP_95, '-o', 'LineWidth', 2); hold on;
ylabel('PICP');
yyaxis right
plot(sens_prob_diff.UncScale, sens_prob_diff.CRPS_Hz, '-s', 'LineWidth', 2);
ylabel('CRPS (Hz)');
grid on;
xlabel('Uncertainty scale');
legend('PICP', 'CRPS', 'Location', 'best');
fontsize(20, "points");


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
plot(convt, (freqsys3(2:end)+1)*60, 'Color', 'k', 'LineWidth', 2); hold on;
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
yline(min((freqsys3(2:end)+1)*60), 'k','LineWidth', 2)
ylabel('Frequency Nadir'); xlabel('# of bins'); fontsize(24,"points");

figure;
h = histogram(min_freq, 30, 'Orientation', 'vertical','Normalization','pdf', 'facecolor', 'r', 'FaceAlpha', 0.5);
hold on
plot(x,y, 'LineWidth', 2)
xline(min((freqsys3(2:end)+1)*60), 'k','LineWidth', 2)
ylabel('Frequency Nadir'); xlabel('# of bins'); fontsize(24,"points");

%% Local Functions

function [x2, y2] = runxy(x, u, A, B, C, D, dt)
   drift = A*x + B*u;
   x2 = x + drift*dt;
   y2 = C*x2 + D*u;
end

function [x2, y2] = runsde(x, u, A, B, C, D, dt,t,unc_scale)
   if nargin < 9
       unc_scale = 1.0;
   end
   drift = A*x + B*u;
   diffusion_f = unc_scale*(0.3*x*(rand()-0.5))*(t>1);
   noise_f = sqrt(dt)*rand(size(drift));
   diffusion_p = unc_scale*(0.15*u*(rand()-0.5))*(t>1);
   noise_p = sqrt(dt)*rand(size(drift));
   x2 = x + drift*dt+ diffusion_p*noise_p;
   y2 = C*x2 + D*u + diffusion_f*noise_f;
end

function proposed = propvoltage(dgenvolts, H_G, u_g, u_l, s, v, unc_scale)
    if nargin < 7
        unc_scale = 1.0;
    end
    sigRow = unc_scale*0.03;        % row scaling
    sigCol = unc_scale*0.03;        % col scaling
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

%% Quantitative Metric Utility Functions Added for Revision

function Tdet = calc_det_metrics(ref_f, pred_f, time_vec, case_name, f_lim)
    ref_f = ref_f(:);
    pred_f = pred_f(:);
    time_vec = time_vec(:);

    L = min([numel(ref_f), numel(pred_f), numel(time_vec)]);
    ref_f = ref_f(1:L);
    pred_f = pred_f(1:L);
    time_vec = time_vec(1:L);

    valid = isfinite(ref_f) & isfinite(pred_f) & isfinite(time_vec);
    ref_f = ref_f(valid);
    pred_f = pred_f(valid);
    time_vec = time_vec(valid);

    err = pred_f - ref_f;

    rmse_f = sqrt(mean(err.^2));
    mae_f = mean(abs(err));
    maxae_f = max(abs(err));

    [ref_nadir, iref] = min(ref_f);
    [pred_nadir, ipred] = min(pred_f);

    nadir_err = abs(pred_nadir - ref_nadir);
    t_nadir_err = abs(time_vec(ipred) - time_vec(iref));

    ss_start = max(1, round(0.8*numel(ref_f)));
    ss_err = mean(pred_f(ss_start:end) - ref_f(ss_start:end));

    vio_ref = ref_nadir < f_lim;
    vio_pred = pred_nadir < f_lim;

    Tdet = table( ...
        string(case_name), rmse_f, mae_f, maxae_f, ...
        ref_nadir, pred_nadir, nadir_err, ...
        time_vec(iref), time_vec(ipred), t_nadir_err, ...
        ss_err, vio_ref, vio_pred, ...
        'VariableNames', { ...
        'Case', 'RMSE_Hz', 'MAE_Hz', 'MaxAE_Hz', ...
        'RefNadir_Hz', 'PredNadir_Hz', 'NadirError_Hz', ...
        'RefTimeToNadir_s', 'PredTimeToNadir_s', 'TimeToNadirError_s', ...
        'SteadyStateMeanError_Hz', 'RefViolation', 'PredViolation'});
end

function Tprob = calc_prob_metrics(ref_f, ens_f, time_vec, case_name, f_lim)
    ref_f = ref_f(:);
    time_vec = time_vec(:);

    if size(ens_f,1) ~= numel(ref_f) && size(ens_f,2) == numel(ref_f)
        ens_f = ens_f.';
    end

    L = min([numel(ref_f), size(ens_f,1), numel(time_vec)]);
    ref_f = ref_f(1:L);
    ens_f = ens_f(1:L,:);
    time_vec = time_vec(1:L);

    valid = isfinite(ref_f) & isfinite(time_vec) & all(isfinite(ens_f),2);
    ref_f = ref_f(valid);
    ens_f = ens_f(valid,:);
    time_vec = time_vec(valid);

    q025 = prctile(ens_f, 2.5, 2);
    q050 = prctile(ens_f, 50, 2);
    q975 = prctile(ens_f, 97.5, 2);

    picp = mean(ref_f >= q025 & ref_f <= q975);
    pinaw = mean(q975 - q025) / (max(ref_f) - min(ref_f) + eps);

    rmse_mean = sqrt(mean((mean(ens_f,2) - ref_f).^2));
    rmse_median = sqrt(mean((q050 - ref_f).^2));
    mae_median = mean(abs(q050 - ref_f));

    crps_t = ensemble_crps_fast(ens_f, ref_f);
    crps_mean = mean(crps_t);

    [ref_nadir, iref] = min(ref_f);
    ens_nadir = min(ens_f, [], 1);

    nadir_mean = mean(ens_nadir);
    nadir_std = std(ens_nadir);
    nadir_q025 = prctile(ens_nadir, 2.5);
    nadir_q05 = prctile(ens_nadir, 5);
    nadir_q50 = prctile(ens_nadir, 50);
    nadir_q95 = prctile(ens_nadir, 95);
    nadir_q975 = prctile(ens_nadir, 97.5);

    nadir_coverage95 = double(ref_nadir >= nadir_q025 && ref_nadir <= nadir_q975);
    nadir_coverage90 = double(ref_nadir >= nadir_q05 && ref_nadir <= nadir_q95);

    violation_rate = mean(ens_nadir < f_lim);
    ref_violation = double(ref_nadir < f_lim);

    interval_width_at_ref_nadir = q975(iref) - q025(iref);

    Tprob = table( ...
        string(case_name), ...
        rmse_mean, rmse_median, mae_median, ...
        picp, pinaw, crps_mean, ...
        ref_nadir, nadir_mean, nadir_std, ...
        nadir_q025, nadir_q05, nadir_q50, nadir_q95, nadir_q975, ...
        nadir_coverage90, nadir_coverage95, ...
        violation_rate, ref_violation, interval_width_at_ref_nadir, ...
        'VariableNames', { ...
        'Case', ...
        'RMSE_Mean_Hz', 'RMSE_Median_Hz', 'MAE_Median_Hz', ...
        'PICP_95', 'PINAW_95', 'CRPS_Hz', ...
        'RefNadir_Hz', 'NadirMean_Hz', 'NadirStd_Hz', ...
        'NadirQ025_Hz', 'NadirQ05_Hz', 'NadirQ50_Hz', 'NadirQ95_Hz', 'NadirQ975_Hz', ...
        'NadirCoverage90', 'NadirCoverage95', ...
        'ViolationRate', 'RefViolation', 'IntervalWidthAtRefNadir_Hz'});
end

function crps_t = ensemble_crps_fast(ens_f, ref_f)
    [T, M] = size(ens_f);
    crps_t = zeros(T,1);

    for k = 1:T
        x = sort(ens_f(k,:));
        y = ref_f(k);

        term1 = mean(abs(x - y));

        i = 1:M;
        coef = 2*i - M - 1;
        term2 = sum(coef .* x) / (M^2);

        crps_t(k) = term1 - term2;
    end
end

function print_metric_table(T, title_str)
    fprintf('\n============================================================\n');
    fprintf('%s\n', title_str);
    fprintf('============================================================\n');
    disp(T);
end

%% Additional Baseline and Voltage-Metric Utility Functions Added for Revision 2

function Tv = calc_voltage_det_metrics(ref_v, pred_v, time_vec, case_name)
    L = min([size(ref_v,1), size(pred_v,1), numel(time_vec)]);
    ref_v = ref_v(1:L,:);
    pred_v = pred_v(1:L,:);
    time_vec = time_vec(1:L);

    E = pred_v - ref_v;
    AE = abs(E);

    mae_v = mean(AE, "all");
    rmse_v = sqrt(mean(E.^2, "all"));
    maxae_v = max(AE, [], "all");

    [ref_min_v, ref_lin_idx] = min(ref_v(:));
    [pred_min_v, pred_lin_idx] = min(pred_v(:));
    [ref_t_idx, ref_bus_idx] = ind2sub(size(ref_v), ref_lin_idx);
    [pred_t_idx, pred_bus_idx] = ind2sub(size(pred_v), pred_lin_idx);

    min_v_error = abs(pred_min_v - ref_min_v);
    time_to_min_v_error = abs(time_vec(pred_t_idx) - time_vec(ref_t_idx));

    post_start = max(1, round(0.8*L));
    post_recovery_mae = mean(AE(post_start:end,:), "all");

    Tv = table( ...
        string(case_name), mae_v, rmse_v, maxae_v, ...
        ref_min_v, pred_min_v, min_v_error, ...
        time_vec(ref_t_idx), time_vec(pred_t_idx), time_to_min_v_error, ...
        ref_bus_idx, pred_bus_idx, post_recovery_mae, ...
        'VariableNames', { ...
        'Case', 'VoltageMAE_pu', 'VoltageRMSE_pu', 'VoltageMaxAE_pu', ...
        'RefMinVoltage_pu', 'PredMinVoltage_pu', 'MinVoltageError_pu', ...
        'RefTimeToMinVoltage_s', 'PredTimeToMinVoltage_s', 'TimeToMinVoltageError_s', ...
        'RefMinVoltageBusIndex', 'PredMinVoltageBusIndex', 'PostRecoveryVoltageMAE_pu'});
end

function lv_conv = convV_from_dgen(dgenvolts, vsens)
    lv_conv = zeros(size(dgenvolts,1), 29);
    for kk = 1:size(dgenvolts,1)
        lv_conv(kk,:) = vsens(:,1:29)' * dgenvolts(kk,:)';
    end
end

function [freq_no_voltage_rec, freq_conv_v_rec, freq_full_sde_rec] = run_probabilistic_baseline_comparison( ...
    Nrsim, N, t, dt, Sbase, cont_size, P_m, ...
    Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
    GFLbase, RGFL, Tconv, KGFL, ...
    ini_load_v, ini_load_p, DER_nominal, v_conv_bus, v_det_bus, unc_scale)

    freq_no_voltage_rec = zeros(N, Nrsim);
    freq_conv_v_rec = zeros(N, Nrsim);
    freq_full_sde_rec = zeros(N, Nrsim);

    % Nominal settings for the Full-SDE baseline. In this baseline, the
    % deterministic SFR drift is fixed and all probabilistic dispersion is
    % represented only by the SDE residual.
    btm_nominal = ini_load_p/sum(ini_load_p) * DER_nominal;
    load_p_nominal = ini_load_p + btm_nominal;
    load_a_nominal = 0.4;
    load_b_nominal = 0.4;
    load_c_nominal = 0.2;
    GFLsize_nominal = GFLbase + DER_nominal;
    GFL_nominal = GFLsize_nominal/Sbase/RGFL * tf(1, [Tconv, 1]) * tf(1, [Tconv, 1]) * tf(KGFL, [1, 0]) / (1 + tf(KGFL, [1, 0]));
    [AGFL_nominal, BGFL_nominal, CGFL_nominal, DGFL_nominal] = ssdata(GFL_nominal);

    for sim = 1:Nrsim
        % Common random uncertainty for conventional baselines except Full-SDE.
        DER_size = DER_nominal * (1 + unc_scale*(rand()-0.5));

        if unc_scale == 0
            btm_load = ini_load_p/sum(ini_load_p) * DER_size;
        else
            btm_load = rand(size(ini_load_v));
            btm_load = btm_load/sum(btm_load) * DER_size;
        end

        load_p = ini_load_p + btm_load;
        load_p = load_p * (1 + unc_scale*(rand()-0.5)*0.2);

        load_a = 0.4 + unc_scale*0.2*(rand()-0.5);
        load_b = 0.4 + unc_scale*0.2*(rand()-0.5);
        load_c = 1 - load_a - load_b;

        GFLsize = GFLbase + DER_size;
        GFL = GFLsize/Sbase/RGFL * tf(1, [Tconv, 1]) * tf(1, [Tconv, 1]) * tf(KGFL, [1, 0]) / (1 + tf(KGFL, [1, 0]));
        [AGFL, BGFL, CGFL, DGFL] = ssdata(GFL);

        % Baseline 1: SFR without voltage recovery or voltage-dependent load relief.
        freq_no_voltage_rec(:,sim) = simulate_sfr_one_baseline( ...
            N, t, dt, Sbase, cont_size, P_m, ...
            Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
            AGFL, BGFL, CGFL, DGFL, GFLsize, ...
            load_p, btm_load, ini_load_v, load_a, load_b, load_c, ...
            [], false, false, true, unc_scale);

        % Baseline 2: SFR using conventional voltage recovery.
        freq_conv_v_rec(:,sim) = simulate_sfr_one_baseline( ...
            N, t, dt, Sbase, cont_size, P_m, ...
            Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
            AGFL, BGFL, CGFL, DGFL, GFLsize, ...
            load_p, btm_load, ini_load_v, load_a, load_b, load_c, ...
            v_conv_bus, true, true, true, unc_scale);

        % Baseline 3: Full-SDE baseline without explicit individual uncertainty modeling.
        freq_full_sde_rec(:,sim) = simulate_sfr_one_baseline( ...
            N, t, dt, Sbase, cont_size, P_m, ...
            Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
            AGFL_nominal, BGFL_nominal, CGFL_nominal, DGFL_nominal, GFLsize_nominal, ...
            load_p_nominal, btm_nominal, ini_load_v, load_a_nominal, load_b_nominal, load_c_nominal, ...
            v_det_bus, true, false, true, unc_scale);
    end
end

function freq = simulate_sfr_one_baseline( ...
    N, t, dt, Sbase, cont_size, P_m, ...
    Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
    AGFL, BGFL, CGFL, DGFL, GFLsize, ...
    load_p, btm_load, ini_load_v, load_a, load_b, load_c, ...
    vmat, use_voltage, use_der_trip, use_sde, unc_scale)

    x_gen = zeros(size(Agen, 1), N);
    x_sim = zeros(size(Anum, 1), N);
    x_gfl = zeros(size(AGFL, 1), N);

    y_gen = zeros(1,N);
    y_gfl = zeros(1,N);
    y_sim = zeros(1,N);
    lc = zeros(1,N);
    trip_size = zeros(1,N);

    for i = 2:N
        y_all = y_gen(:, i-1) + y_gfl(:, i-1) + lc(:, i-1);
        % u = P_m(t(i-1)) - y_all;
        p_dist = (-cont_size/Sbase) * (t(i-1) >= 1);
        u = p_dist - y_all;

        if use_sde
            [x_sim(:, i), y_sim(:, i)] = runsde(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt, t(i), unc_scale);
        else
            [x_sim(:, i), y_sim(:, i)] = runxy(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt);
        end

        if use_voltage
            if i <= size(vmat,1)
                v_now = vmat(i,:);
            else
                v_now = vmat(end,:);
            end

            if use_der_trip
                trip_signal = DER_vrt(size(vmat,2), t(i), dt, v_now);
                % trip_size(i) = sum(btm_load .* (trip_signal > 0));
                trip_size(i) = sum(btm_load(:) .* (trip_signal(:) > 0));
            end

            lc(:,i) = sum(load_p .* (load_a.*(v_now./ini_load_v).^2 + load_b.*(v_now./ini_load_v) + load_c - 1)) / (sum(load_p) + cont_size);
        end

        [x_gen(:, i), y_gen(:, i)] = runxy(x_gen(:, i-1), y_sim(:,i), Agen, Bgen, Cgen, Dgen, dt);
        [x_gfl(:, i), y_gfl(:, i)] = runxy(x_gfl(:, i-1), y_sim(:,i), AGFL, BGFL, CGFL, DGFL, dt);
        y_gfl(:, i) = max(-0.4*GFLsize/Sbase, y_gfl(:,i));

        if use_der_trip && GFLsize > 0
            y_gfl(:, i) = y_gfl(:,i) * max((GFLsize - trip_size(i))/GFLsize, 0);
        end
    end

    freq = (1 + y_sim(:)) * 60;
end

function safe_writetable(T, filename)
    try
        writetable(T, filename);
    catch ME
        warning('Could not write %s: %s', filename, ME.message);
    end
end