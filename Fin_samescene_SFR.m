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

%% Probabilistic SFR
Nrsim = 500;

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
        DER_size = 600*(1+(rand()-0.5));        
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

    % Voltage generation
    pro_prop = propvoltage(dgenvolts, H_G, u_g, u_l, s, v) + Volt(1,1:29);
    
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
end

[mean_freq, ci_z, ci_m] = find_ci(freq_rec);
[mean_gen, ci_g, ci_g2] = find_ci(y_gen_rec);
[mean_gfl, ci_i, ci_i2] = find_ci(y_gfl_rec);
[mean_load, ci_l, ci_l2] = find_ci(y_load_rec);
[mean_v, ci_v, ci_v2] = find_ci(rec_prop);

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

function [x2, y2] = runsde(x, u, A, B, C, D, dt,t)
   drift = A*x + B*u;
   diffusion_f = (0.3*x*(rand()-0.5))*(t>1);
   noise_f = sqrt(dt)*rand(size(drift));
   diffusion_p = (0.15*u*(rand()-0.5))*(t>1);
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