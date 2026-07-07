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
Start = 950;
Stop = 1300;
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
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Ud(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(Start:Stop),output_rec(Start:Stop) + D_P_DC, ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Ud(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('D轴开环正向方波响应：\n')
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('加性误差修正前相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm(output_rec(Start:Stop) - output_sim(Start:Stop)) / ...
    norm(output_rec(Start:Stop) - mean(output_rec(Start:Stop))));
fprintf('加性误差修正前拟合优度：%.5f%%\n',r);
cov_sim_rec = cov((output_rec(Start:Stop)+D_P_DC), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var((output_rec(Start:Stop)+D_P_DC)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('加性误差修正后相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm((output_rec(Start:Stop)+D_P_DC) - output_sim(Start:Stop)) / ...
    norm((output_rec(Start:Stop)+D_P_DC) - mean(output_rec(Start:Stop))));
fprintf('加性误差修正后拟合优度：%.5f%%\n\n',r);

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
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Ud(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(Start:Stop),output_rec(Start:Stop) + D_N_DC, ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Ud(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('D轴开环负向方波响应：\n');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('加性误差修正前相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm(output_rec(Start:Stop) - output_sim(Start:Stop)) / ...
    norm(output_rec(Start:Stop) - mean(output_rec(Start:Stop))));
fprintf('加性误差修正前拟合优度：%.5f%%\n',r);
cov_sim_rec = cov((output_rec(Start:Stop)+D_N_DC), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var((output_rec(Start:Stop)+D_N_DC)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('加性误差修正后相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm((output_rec(Start:Stop)+D_N_DC) - output_sim(Start:Stop)) / ...
    norm((output_rec(Start:Stop)+D_N_DC) - mean(output_rec(Start:Stop))));
fprintf('加性误差修正后拟合优度：%.5f%%\n\n',r);


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
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Uq(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(Start:Stop),output_rec(Start:Stop).*Q_P_GAIN, ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Uq(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('Q轴开环正向方波响应：\n');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('乘性误差修正前相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm(output_rec(Start:Stop) - output_sim(Start:Stop)) / ...
    norm(output_rec(Start:Stop) - mean(output_rec(Start:Stop))));
fprintf('乘性误差修正前拟合优度：%.5f%%\n',r);
cov_sim_rec = cov(output_rec(Start:Stop)*Q_P_GAIN, output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var((output_rec(Start:Stop)*Q_P_GAIN)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('乘性误差修正后相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm((output_rec(Start:Stop)*Q_P_GAIN) - output_sim(Start:Stop)) / ...
    norm((output_rec(Start:Stop)*Q_P_GAIN) - mean(output_sim(Start:Stop))));
fprintf('乘性误差修正后拟合优度：%.5f%%\n\n',r);

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
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Uq(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正前');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
subplot(1,2,2);
plot(t(Start:Stop),output_rec(Start:Stop).*Q_N_GAIN, ...
    t(Start:Stop),output_sim(Start:Stop), ...
    t(Start:Stop),Uq(Start:Stop).*0.1,LineWidth=2);
subtitle('误差修正后');
legend('实测输出','仿真输出','输入信号（衰减0.1）');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
fprintf('Q轴开环正向方波响应：\n');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('乘性误差修正前相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm(output_rec(Start:Stop) - output_sim(Start:Stop)) / ...
    norm(output_rec(Start:Stop) - mean(output_rec(Start:Stop))));
fprintf('乘性误差修正前拟合优度：%.5f%%\n',r);
cov_sim_rec = cov(output_rec(Start:Stop)*Q_N_GAIN, output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var((output_rec(Start:Stop)*Q_N_GAIN)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('乘性误差修正后相关系数：%.5f%%\n',r*100);
r = 100 * (1 - norm((output_rec(Start:Stop)*Q_N_GAIN) - output_sim(Start:Stop)) / ...
    norm((output_rec(Start:Stop)*Q_N_GAIN) - mean(output_sim(Start:Stop))));
fprintf('乘性误差修正后拟合优度：%.5f%%\n\n',r);