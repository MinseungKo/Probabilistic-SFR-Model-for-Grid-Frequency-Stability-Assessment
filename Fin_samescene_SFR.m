clear all
close all
clc

%% Import Data
rng('shuffle')

PFile='loadbus_3_500.xlsx';
vsens = csvread('volt_sensitivity.csv'); %% Voltage sensitivity obtained by reactive power perturbation
H = csvread('voltH.csv');

data = xlsread(PFile);
data = data(5:end, :);

time = data(:, 1);
Ps = data(:, 2:11);
Freq = data(:, 12:50);
Volt = data(:, 51:end);

dt = time(2,1);
T = time(end,1);
N = round(T/dt)+1; % Total solve steps

%% Generator & System Informations
Basepow = 100;
Sbase = sum(Ps(1,:))*Basepow; %MW
SG_num = [2, 3, 4, 7, 9, 10];
Genbase = sum(Ps(1,SG_num))*Basepow;
GFLbase = Sbase-Genbase;

Hs = [2.3, 4.9, 4.3, 5.2, 3.15, 5.0]; % 31, 32, 33, 36, 38, 39
freqsys = sum(Hs.*Freq(:, SG_num),2)/sum(Hs);

Hsys = 0;
for i=1:6
    Hsys = Hsys+Hs(i)*Ps(1,SG_num(i))*100;
end
Dsys=0;
Hsys = Hsys/Sbase;

%% Data for Voltage Estimation
ini_load_v = [1.0345, 1.0221, 1.0102, 1.0098, 1.0109, 1.024, 1.0374, ...
                1.0355, 0.9941, 1.0359, 1.0474,1.0434,1.0586,1.0547,1.0428,1.0515,1.0518,0.982,1.03];
ini_load_p = [322, 400, 233.8, 522, 7.5, 220, 329.4, 158, 528, 274, 274.5, 208.6, 224, 139, 181, ...
                206,183.5, 9.2, 1104];

gen_bus = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39];
gen_load_bus = [2, 10];
load_bus = [3,4,7,8,12,15,16,18,20,21,23,24,25,26,27,28,29];

genvolts = Volt(:, gen_bus);
dgenvolts = genvolts-genvolts(1,:);

conventional = convV(Volt, vsens);

% Proposed Method
H_L = H(1:29, :);
H_G = H(30:end, :);

[u, s, v] = svd(H);
M = 10;
u_l = u(1:29, 1:M);
u_g = u(30:end, 1:M);
s = s(1:M, :);

for i=1:size(genvolts,1)
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
lv_conv = conventional + Volt(1,1:29);
lv_prop = ori_prop + Volt(1,1:29);
lv_real = Volt(:, 1:29);

%% Voltage Estimation Figures
Err_prop = abs(lv_real-lv_prop);
Err_conv = abs(lv_real-lv_conv);
mae_prop = mean(mean(Err_prop));
mae_conv = mean(mean(Err_conv));

minlvolts = min(lv_real);
[minlv, minlvind] = mink(minlvolts, 3);
[maxlv, maxlvind] = maxk(minlvolts, 3);
lvind = [minlvind, maxlvind];

cmin = 0;
cmax = max([Err_conv(:); Err_prop(:)]);

figure;
subplot(2, 1,1)
imagesc(time, 1:size(Err_prop, 2), Err_conv')
set(gca,'YDir','normal')
caxis([cmin cmax])
title('Conventional method'); ylabel('Load bus #'); fontsize(24,"points");
colormap(jet)
colorbar

subplot(2, 1,2)
imagesc(time, 1:size(Err_prop, 2), Err_prop')
set(gca,'YDir','normal')
caxis([cmin cmax])
title('Proposed method'); xlabel('Time (s)'); ylabel('Load bus #'); fontsize(24,"points");
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
v_load_bus_conv = lv_conv(:, load_bus);
v_corr_gen_load_bus = genvolts(:, gen_load_bus);
v_load_bus = [v_corr_only_load_bus, v_corr_gen_load_bus];
v_load_bus_comp = [v_load_bus_conv, v_corr_gen_load_bus];

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
    % [x_sim(:, i), y_sim(:, i)] = runsde(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt, t(i));
    % y_sim(:, i) = y_sim(:,i) + Jump(t(i)) + Jump_rec(t(i));

    % volts
    lc(:,i) = load_riv(i, v_load_bus); %volts);%, volts(i,:)); %loadcorrvolts
    % lc(:,i) = load_change(t(i), y_sim(i)*60);
    [x_gen(:, i), y_gen(:, i)] = runxy(x_gen(:, i-1), y_sim(:,i), Agen, Bgen, Cgen, Dgen, dt);
    [x_gfl(:, i), y_gfl(:, i)] = runxy(x_gfl(:, i-1), y_sim(:,i), AGFL, BGFL, CGFL, DGFL, dt);
    y_gfl(:, i) = max(-0.3*GFLsize/Sbase, y_gfl(:,i));
end

x_gen2 = zeros(size(Agen, 1), N);
x_sim2 = zeros(size(Anum, 1), N);
x_gfl2 = zeros(size(AGFL, 1), N);

y_gen2 = zeros(1,N);
y_gfl2 = zeros(1,N);
y_sim2 = zeros(1,N);

lc2 = zeros(1,N);
for i = 2:N   
    y_all2 = y_gen2(:, i-1)+y_gfl2(:, i-1)+lc2(:, i-1);
    u2 = P_m(t(i-1))-y_all2;

    [x_sim2(:, i), y_sim2(:, i)] = runxy(x_sim2(:, i-1), u2, Anum, Bnum, Cnum, Dnum, dt);

    % volts
    lc2(:,i) = load_riv(i, v_load_bus_comp);
    [x_gen2(:, i), y_gen2(:, i)] = runxy(x_gen2(:, i-1), y_sim2(:,i), Agen, Bgen, Cgen, Dgen, dt);
    [x_gfl2(:, i), y_gfl2(:, i)] = runxy(x_gfl2(:, i-1), y_sim2(:,i), AGFL, BGFL, CGFL, DGFL, dt);
    y_gfl2(:, i) = max(-0.3*GFLsize/Sbase, y_gfl2(:,i));
end

x_gen3 = zeros(size(Agen, 1), N);
x_sim3 = zeros(size(Anum, 1), N);
x_gfl3 = zeros(size(AGFL, 1), N);

y_gen3 = zeros(1,N);
y_gfl3 = zeros(1,N);
y_sim3 = zeros(1,N);
for i = 2:N   
    y_all3 = y_gen3(:, i-1)+y_gfl3(:, i-1);
    u3 = P_m(t(i-1))-y_all3;

    [x_sim3(:, i), y_sim3(:, i)] = runxy(x_sim3(:, i-1), u3, Anum, Bnum, Cnum, Dnum, dt);

    % volts
    [x_gen3(:, i), y_gen3(:, i)] = runxy(x_gen3(:, i-1), y_sim3(:,i), Agen, Bgen, Cgen, Dgen, dt);
    [x_gfl3(:, i), y_gfl3(:, i)] = runxy(x_gfl3(:, i-1), y_sim3(:,i), AGFL, BGFL, CGFL, DGFL, dt);
    y_gfl3(:, i) = max(-0.3*GFLsize/Sbase, y_gfl3(:,i));
end
freq_det = (1+y_sim)*60;
freq2 = (1+y_sim2)*60;
freq3 = (1+y_sim3)*60;

%% Probabilistic SFR and Uncertainty Sensitivity
Nrsim = 500;
unc_scale_list = [0, 0.5, 1.0, 1.5, 2.0];
base_unc_scale = 1.0;
sens_prob_same = table();
base_prob_store = struct();

for iu = 1:numel(unc_scale_list)
    unc_scale = unc_scale_list(iu);
    fprintf('Running same-scenario probabilistic SFR with uncertainty scale = %.2f\n', unc_scale);

rec_prop = zeros([size(data,1), 29, Nrsim]);

unc_der_s = 1;
unc_lmp_d = 1;
unc_v_d = 1;
unc_la_s = 1;

freq_rec = zeros(size(data,1)-1, Nrsim);
y_gen_rec = zeros(size(data,1)-1, Nrsim);
y_load_rec = zeros(size(data,1)-1, Nrsim);
y_gfl_rec = zeros(size(data,1)-1, Nrsim);

for sim=1:Nrsim

    % BTM & Pure Load
    if unc_der_s == 1
        DER_size = 600*(1 + unc_scale*(rand()-0.5));        
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

    % Voltage generation
    pro_prop = propvoltage(dgenvolts, H_G, u_g, u_l, s, v, unc_scale) + Volt(1,1:29);

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
end



    ref_freq_sens = (freqsys(2:end) + 1) * 60;
    time_sens = t(:);
    Lsens = min([numel(ref_freq_sens), size(freq_rec,1), numel(time_sens)]);
    Tprob_i = calc_prob_metrics(ref_freq_sens(1:Lsens), freq_rec(1:Lsens,:), time_sens(1:Lsens), "Same-Uncertainty sensitivity", 59.5);
    Tprob_i.UncScale = unc_scale;
    Tprob_i = movevars(Tprob_i, 'UncScale', 'After', 'Case');
    sens_prob_same = [sens_prob_same; Tprob_i];

    if abs(unc_scale - base_unc_scale) < 1e-12
        base_prob_store.freq_rec = freq_rec;
        base_prob_store.y_gen_rec = y_gen_rec;
        base_prob_store.y_gfl_rec = y_gfl_rec;
        base_prob_store.y_load_rec = y_load_rec;
        base_prob_store.rec_prop = rec_prop;
        base_prob_store.min_freq = min_freq;
    end
end

freq_rec = base_prob_store.freq_rec;
y_gen_rec = base_prob_store.y_gen_rec;
y_gfl_rec = base_prob_store.y_gfl_rec;
y_load_rec = base_prob_store.y_load_rec;
rec_prop = base_prob_store.rec_prop;
min_freq = base_prob_store.min_freq;



%% Probabilistic Baseline Comparison: Same Scenario
[freq_rec_novolt, freq_rec_conv, freq_rec_fullsde] = run_probabilistic_baseline_comparison( ...
    Nrsim, N, t, dt, Sbase, cont_size, P_m, ...
    Anum, Bnum, Cnum, Dnum, Agen, Bgen, Cgen, Dgen, ...
    GFLbase, RGFL, Tconv, KGFL, ...
    ini_load_v, ini_load_p, 600, v_load_bus_comp, v_load_bus, base_unc_scale);

[mean_freq, ci_z, ci_m] = find_ci(freq_rec);
[mean_gen, ci_g, ci_g2] = find_ci(y_gen_rec);
[mean_gfl, ci_i, ci_i2] = find_ci(y_gfl_rec);
[mean_load, ci_l, ci_l2] = find_ci(y_load_rec);
[mean_v, ci_v, ci_v2] = find_ci(rec_prop);


%% Quantitative Metrics: Same Scenario
f_lim = 59.5;

ref_freq_same = (freqsys(2:end) + 1) * 60;
time_metric_same = t(:);

Lmet = min([numel(ref_freq_same), numel(freq_det), numel(freq2), numel(freq3), size(freq_rec,1), numel(time_metric_same)]);
ref_freq_same = ref_freq_same(1:Lmet);
time_metric_same = time_metric_same(1:Lmet);

freq_det_m = freq_det(1:Lmet);
freq_conv_m = freq2(1:Lmet);
freq_novolt_m = freq3(1:Lmet);
freq_rec_m = freq_rec(1:Lmet,:);

Tdet_same = [
    calc_det_metrics(ref_freq_same, freq_det_m, time_metric_same, "Same-Proposed deterministic SFR", f_lim);
    calc_det_metrics(ref_freq_same, freq_conv_m, time_metric_same, "Same-Conventional voltage recovery SFR", f_lim);
    calc_det_metrics(ref_freq_same, freq_novolt_m, time_metric_same, "Same-No voltage-dependent load relief", f_lim)
];

Vref_same = lv_real;
Vprop_same = lv_prop;
Vconv_same = lv_conv;

Lv = min([size(Vref_same,1), size(Vprop_same,1), size(Vconv_same,1)]);
Vref_same = Vref_same(1:Lv,:);
Vprop_same = Vprop_same(1:Lv,:);
Vconv_same = Vconv_same(1:Lv,:);

voltage_metrics_same = table( ...
    ["Conventional voltage recovery"; "Proposed voltage recovery"], ...
    [mean(abs(Vconv_same - Vref_same), "all"); mean(abs(Vprop_same - Vref_same), "all")], ...
    [sqrt(mean((Vconv_same - Vref_same).^2, "all")); sqrt(mean((Vprop_same - Vref_same).^2, "all"))], ...
    [max(abs(Vconv_same - Vref_same), [], "all"); max(abs(Vprop_same - Vref_same), [], "all")], ...
    'VariableNames', {'Case', 'VoltageMAE_pu', 'VoltageRMSE_pu', 'VoltageMaxAE_pu'} ...
);

voltage_det_metrics_same = [
    calc_voltage_det_metrics(Vref_same, Vconv_same, time(1:Lv), "Same-Conventional voltage recovery");
    calc_voltage_det_metrics(Vref_same, Vprop_same, time(1:Lv), "Same-Proposed voltage recovery")
];

Tprob_same = calc_prob_metrics(ref_freq_same, freq_rec_m, time_metric_same, "Same-Proposed probabilistic SFR", f_lim);


freq_rec_novolt_m = freq_rec_novolt(1:Lmet,:);
freq_rec_conv_m = freq_rec_conv(1:Lmet,:);
freq_rec_fullsde_m = freq_rec_fullsde(1:Lmet,:);

Tprob_compare_same = [
    Tprob_same;
    calc_prob_metrics(ref_freq_same, freq_rec_novolt_m, time_metric_same, "Same-Baseline SFR without voltage recovery", f_lim);
    calc_prob_metrics(ref_freq_same, freq_rec_conv_m, time_metric_same, "Same-Baseline SFR with conventional voltage recovery", f_lim);
    calc_prob_metrics(ref_freq_same, freq_rec_fullsde_m, time_metric_same, "Same-Baseline full-SDE SFR without explicit uncertainty modeling", f_lim)
];

rec_prop_m = rec_prop(1:min(size(rec_prop,1), size(lv_real,1)),:,:);
Vref_prob = lv_real(1:size(rec_prop_m,1),:);

qV025 = prctile(rec_prop_m, 2.5, 3);
qV975 = prctile(rec_prop_m, 97.5, 3);

V_PICP_95 = mean((Vref_prob >= qV025) & (Vref_prob <= qV975), "all");
V_PINAW_95 = mean(qV975 - qV025, "all") / (max(Vref_prob, [], "all") - min(Vref_prob, [], "all") + eps);

voltage_prob_metrics_same = table( ...
    "Same-Probabilistic voltage recovery", V_PICP_95, V_PINAW_95, ...
    'VariableNames', {'Case', 'VoltagePICP_95', 'VoltagePINAW_95'} ...
);

print_metric_table(Tdet_same, "Deterministic frequency metrics: same scenario");
print_metric_table(voltage_metrics_same, "Deterministic voltage recovery metrics: same scenario");
print_metric_table(voltage_det_metrics_same, "Expanded deterministic voltage recovery metrics: same scenario");
print_metric_table(Tprob_same, "Probabilistic frequency metrics: same scenario");
print_metric_table(Tprob_compare_same, "Probabilistic frequency metric comparison: same scenario");
print_metric_table(voltage_prob_metrics_same, "Probabilistic voltage coverage metrics: same scenario");
print_metric_table(sens_prob_same, "Uncertainty sensitivity metrics: same scenario");

safe_writetable(Tdet_same, "metrics_same_deterministic_frequency.csv");
safe_writetable(voltage_metrics_same, "metrics_same_voltage.csv");
safe_writetable(voltage_det_metrics_same, "metrics_same_voltage_expanded.csv");
safe_writetable(Tprob_same, "metrics_same_probabilistic_frequency.csv");
safe_writetable(Tprob_compare_same, "metrics_same_probabilistic_frequency_comparison.csv");
safe_writetable(voltage_prob_metrics_same, "metrics_same_probabilistic_voltage.csv");
safe_writetable(sens_prob_same, "metrics_same_uncertainty_sensitivity.csv");

figure;
plot(sens_prob_same.UncScale, sens_prob_same.NadirQ05_Hz, '-o', 'LineWidth', 2); hold on;
plot(sens_prob_same.UncScale, sens_prob_same.NadirQ50_Hz, '-s', 'LineWidth', 2);
plot(sens_prob_same.UncScale, sens_prob_same.NadirQ95_Hz, '-^', 'LineWidth', 2);
yline(f_lim, 'k--', 'LineWidth', 2);
grid on;
xlabel('Uncertainty scale'); ylabel('Frequency nadir (Hz)');
legend('5% nadir', 'Median nadir', '95% nadir', 'Security limit', 'Location', 'best');
fontsize(20, "points");

figure;
yyaxis left
plot(sens_prob_same.UncScale, sens_prob_same.PICP_95, '-o', 'LineWidth', 2); hold on;
ylabel('PICP');
yyaxis right
plot(sens_prob_same.UncScale, sens_prob_same.CRPS_Hz, '-s', 'LineWidth', 2);
ylabel('CRPS (Hz)');
grid on;
xlabel('Uncertainty scale');
legend('PICP', 'CRPS', 'Location', 'best');
fontsize(20, "points");


convt = t';
figure;
fill([convt; flipud(convt)], [ci_m(:,1); flipud(ci_m(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'r', 'EdgeAlpha', 0.15); hold on;
fill([convt; flipud(convt)], [ci_z(:,1); flipud(ci_z(:,2))], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.3);
plot(convt, (freqsys(2:end)+1)*60, 'Color', 'k', 'LineWidth', 2); hold on;
plot(convt, freq_det, 'LineStyle', ':', 'LineWidth', 2, 'Color', 'r');
plot(t, freq2, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-.');
plot(t, freq3, 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--');
grid on; xlabel('Time (s)'); ylabel('Frequency (Hz)'); fontsize(24,"points");


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

figure(7)
h = histogram(min_freq, 30, 'Orientation', 'horizontal','Normalization','pdf', 'facecolor', 'r', 'FaceAlpha', 0.5);
hold on
plot(y,x, 'r','LineWidth', 2)
yline(min((freqsys(2:end)+1)*60), 'k','LineWidth', 2)
ylabel('Frequency Nadir'); xlabel('# of bins'); fontsize(24,"points");

figure(8)
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

function conventional = convV(Volt, vsens)
    allbusnum = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39];
    
    genvolts = Volt(:, allbusnum);
    dgenvolts = genvolts-genvolts(1,:);
    
    % Conventional Method
    for i=1:size(genvolts,1)
        conventional(i,:) = vsens(:, 1:29)'*dgenvolts(i,:)';
    end
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

