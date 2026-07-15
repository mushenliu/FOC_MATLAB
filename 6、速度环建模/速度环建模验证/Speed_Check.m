close all;
clear;
clc;
s = tf('s');
%% id=0情况下的理论模型
Rs = 5.70125;
Ld = 1.86152e-3;
Psi = 9.65548e-3;
Lq = 2.16828e-3;
p = 14;
f = 1000;
Ts = 1/f;
%速度环零极点和开环增益
p1 = -606.19068;
p2 = -11.11839;
z1 = 1861.76409;
z2 = -2629.38827;
K_m = 15066.81290;

G_RL = 1/(Rs + Lq * s);
K = 1.5 * p * Psi;
G_m = (s-z1)*(s-z2)/((s-p1)*(s-p2));
G_m = G_m * K_m / dcgain(G_m);
Kf = pi / 30 * p * Psi;
G_C = feedback(G_RL*K*G_m,Kf);
% G_C = c2d(G_C,Ts,'zoh');

%阶跃统计序列段
Start = 950;
Stop = 1150;

%% 正向方波
data =  readmatrix("Speed_Step_P.txt", 'NumHeaderLines', 1);
Uq = data(end-2046:end,9);
Ud = data(end-2046:end,8);
theta_e = data(end-2046:end,2);
Id = data(end-2046:end,5);
Iq = data(end-2046:end,6);
n = data(end-2046:end,3);
we = n / 60 * 2 * pi * p;
Teqd = 1.5*p.*Iq.*((Ld-Lq).*Id+Psi);
Teq = 1.5*p*Psi.*Iq;
Eqd = we.*(Psi+Ld.*Id);
Eq = we.*Psi;
figure(Name='理想假设分析');
subplot(1,2,1);
plot(Teq);
hold on;
plot(Teqd);
grid on;
grid minor;
subtitle('电磁转矩分析');
ylabel('电磁转矩/N·m');
legend('id=0电磁转矩','id≠0电磁转矩');
subplot(1,2,2);
plot(Eq);
hold on;
plot(Eqd);
grid on;
grid minor;
ylabel('反电动势/V');
subtitle('反电动势分析');
legend('id=0反电动势','id≠0反电动势');
fprintf('电磁转矩理想化分析：\n');
Signal_Analyse(Teq,Teqd);
fprintf('反电动势分析：\n');
Signal_Analyse(Eq,Eqd);
t = Ts:Ts:Ts*(length(Uq));
output_sim = lsim(G_C,Uq(1:end),t);
figure(Name='正向');
plot(t(Start:Stop),n(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Uq(Start:Stop)*10,LineWidth=2);
legend('实测输出','仿真输出','输入信号（增益10倍）');
grid on;
grid minor;
ylabel('转速/rpm');
xlabel('时间/s');
fprintf('正向方波响应：\n');
Signal_Analyse(n(Start:Stop),output_sim(Start:Stop));

%% 负向方波
data =  readmatrix("Speed_Step_N.txt", 'NumHeaderLines', 1);
Uq = data(end-2046:end,9);
Ud = data(end-2046:end,8);
theta_e = data(end-2046:end,2);
Id = data(end-2046:end,5);
Iq = data(end-2046:end,6);
n = data(end-2046:end,3);
we = n / 60 * 2 * pi * p;
Teqd = 1.5*p.*Iq.*((Ld-Lq).*Id+Psi);
Teq = 1.5*p*Psi.*Iq;
Eqd = we.*(Psi+Ld.*Id);
Eq = we.*Psi;
figure(Name='理想假设分析');
subplot(1,2,1);
plot(Teq);
hold on;
plot(Teqd);
grid on;
grid minor;
subtitle('电磁转矩分析');
legend('id=0电磁转矩','id≠0电磁转矩');
subplot(1,2,2);
plot(Eq);
grid on;
grid minor;
hold on;
plot(Eqd);
subtitle('反电动势分析');
legend('id=0反电动势','id≠0反电动势');
fprintf('\n电磁转矩理想化分析：\n');
Signal_Analyse(Teq,Teqd);
fprintf('反电动势分析：\n');
Signal_Analyse(Eq,Eqd);
output_sim = lsim(G_C,Uq(1:end),t);
figure(Name='负向');
plot(t(Start:Stop),n(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Uq(Start:Stop)*10,LineWidth=2);
legend('实测输出','仿真输出','输入信号（增益10倍）');
grid on;
grid minor;
ylabel('转速/rpm');
xlabel('时间/s');
fprintf('负向方波响应：\n')
Signal_Analyse(n(Start:Stop),output_sim(Start:Stop));