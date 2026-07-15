close all;
clear;
clc;
s=tf('s');
%% 仿真模型构建
Rs = 5.70125;
Ld = 1.86152e-3;
Psi = 9.65548e-3;
Lq = 2.16828e-3;
Ts = 1/(170000000 / 5665 / 2);
%滤波器系数
fi = 0.5;
%偏置修正（D轴加性误差，Q轴乘性误差）
D_P_DC = -0.02;
D_N_DC = -0.015;
Q_P_GAIN = 1.1;
Q_N_GAIN = 1.1;

G_D = c2d(1/(Rs+Ld*s) * exp(-1*Ts*s) * exp(-260e-9 * s) * 1/(120e-9 * s + 1), ...
    Ts,'zoh');
G_Q = c2d(1/(Rs+Lq*s) * exp(-1*Ts*s) * exp(-260e-9 * s) * 1/(120e-9 * s + 1), ...
    Ts,'zoh');

%阶跃统计序列段
D_Start = 950;
D_Stop = 1300;
Q_Start = 950;
Q_Stop = 1300;
%% D轴正向方波
data =  readmatrix("D_Step_P.txt", 'NumHeaderLines', 1);
Ud = data(end-2046:end,8);
theta_e = data(end-2046:end,2);
Id = data(end-2046:end,5);
we = zeros(length(Id),1);
we(1) = 0;
for i=2:length(we)
    we(i)=(theta_e(i)-theta_e(i-1))/Ts * 2 * pi / 360;
    if(abs(theta_e(i)-theta_e(i-1))>50)
        we(i)=we(i-1);
    end
end
we_lp = lowpass(we,fi);
input = Ud + we_lp .* (Lq .* data(end-2046:end,6));
output_rec = data(end-2046:end,5);
t = Ts:Ts:Ts*(length(input));
output_sim = lsim(G_D,input(1:end),t);
figure(Name='D轴正向');
subplot(1,2,1);
plot(t(D_Start:D_Stop),output_rec(D_Start:D_Stop), ...
    t(D_Start:D_Stop),output_sim(D_Start:D_Stop), ...
    t(D_Start:D_Stop),Ud(D_Start:D_Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(D_Start:D_Stop),output_rec(D_Start:D_Stop) + D_P_DC, ...
    t(D_Start:D_Stop),output_sim(D_Start:D_Stop), ...
    t(D_Start:D_Stop),Ud(D_Start:D_Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('D轴开环正向方波响应（加性误差修正前）：\n')
Signal_Analyse(output_rec(D_Start:D_Stop),output_sim(D_Start:D_Stop));
fprintf('D轴开环正向方波响应（加性误差修正后）：\n')
Signal_Analyse(output_rec(D_Start:D_Stop)+D_P_DC,output_sim(D_Start:D_Stop));


%% D轴负向方波
data =  readmatrix("D_Step_N.txt", 'NumHeaderLines', 1);
Ud = data(end-2046:end,8);
theta_e = data(end-2046:end,2);
Id = data(end-2046:end,5);
we = zeros(length(Id),1);
we(1) = 0;
for i=2:length(we)
    we(i)=(theta_e(i)-theta_e(i-1))/Ts * 2 * pi / 360;
    if(abs(theta_e(i)-theta_e(i-1))>50)
        we(i)=we(i-1);
    end
end
we_lp = lowpass(we,fi);
input = Ud + we_lp .* (Lq .* data(end-2046:end,6));
output_rec = data(end-2046:end,5);
output_sim = lsim(G_D,input(1:end),t);
figure(Name='D轴负向');
subplot(1,2,1);
plot(t(D_Start:D_Stop),output_rec(D_Start:D_Stop), ...
    t(D_Start:D_Stop),output_sim(D_Start:D_Stop), ...
    t(D_Start:D_Stop),Ud(D_Start:D_Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(D_Start:D_Stop),output_rec(D_Start:D_Stop) + D_N_DC, ...
    t(D_Start:D_Stop),output_sim(D_Start:D_Stop), ...
    t(D_Start:D_Stop),Ud(D_Start:D_Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nD轴开环负向方波响应（加性误差修正前）：\n')
Signal_Analyse(output_rec(D_Start:D_Stop),output_sim(D_Start:D_Stop));
fprintf('D轴开环负向方波响应（加性误差修正后）：\n')
Signal_Analyse(output_rec(D_Start:D_Stop)+D_P_DC,output_sim(D_Start:D_Stop));



%% Q轴正向方波
data = readmatrix("Q_Step_P.txt", 'NumHeaderLines', 1);
Uq = data(end-2046:end,9);
theta_e = data(end-2046:end,2);
Id = data(end-2046:end,5);
we = zeros(length(Id),1);
we(1) = 0;
for i=2:length(we)
    we(i)=(theta_e(i)-theta_e(i-1))/Ts * 2 * pi / 360;
    if(abs(theta_e(i)-theta_e(i-1))>50)
        we(i)=we(i-1);
    end
end
we_lp = lowpass(we,fi);
input = Uq - we_lp .* (Ld * Id + Psi);
output_rec = data(end-2046:end,6);
output_sim = lsim(G_Q,input(1:end),t);
figure(Name='Q轴正向');
subplot(1,2,1);
plot(t(Q_Start:Q_Stop),output_rec(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),Uq(Q_Start:Q_Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(Q_Start:Q_Stop),output_rec(Q_Start:Q_Stop).*Q_P_GAIN, ...
    t(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),Uq(Q_Start:Q_Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nQ轴开环正向方波响应（乘性误差修正前）：\n')
Signal_Analyse(output_rec(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop));
fprintf('Q轴开环正向方波响应（乘性误差修正后）：\n')
Signal_Analyse(output_rec(Q_Start:Q_Stop).*Q_P_GAIN,output_sim(Q_Start:Q_Stop));


%% Q轴负向阶跃
data = readmatrix("Q_Step_N.txt", 'NumHeaderLines', 1);
Uq = data(end-2046:end,9);
theta_e = data(end-2046:end,2);
Id = data(end-2046:end,5);
we = zeros(length(Id),1);
we(1) = 0;
for i=2:length(we)
    we(i)=(theta_e(i)-theta_e(i-1))/Ts * 2 * pi / 360;
    if(abs(theta_e(i)-theta_e(i-1))>50)
        we(i)=we(i-1);
    end
end
we_lp = lowpass(we,fi);
input = Uq - we_lp .* (Ld * Id + Psi);
output_rec = data(end-2046:end,6);
output_sim = lsim(G_Q,input(1:end),t);
figure(Name='Q轴负向');
subplot(1,2,1);
plot(t(Q_Start:Q_Stop),output_rec(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),Uq(Q_Start:Q_Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(Q_Start:Q_Stop),output_rec(Q_Start:Q_Stop).*Q_N_GAIN, ...
    t(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),Uq(Q_Start:Q_Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nQ轴开环负向方波响应（乘性误差修正前）：\n')
Signal_Analyse(output_rec(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop));
fprintf('Q轴开环负向方波响应（乘性误差修正后）：\n')
Signal_Analyse(output_rec(Q_Start:Q_Stop).*Q_P_GAIN,output_sim(Q_Start:Q_Stop));