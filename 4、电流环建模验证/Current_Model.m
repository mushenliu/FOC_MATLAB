close all;
clear;
clc;
s=tf('s');
%% 仿真模型构建
Rs = 5.70125;
Ld = 0.00184357;
Lq = 0.00194861;
Ts = 1/(170000000 / 5665 / 2);
Psi = 0.009600568;
%滤波器系数
fi = 0.001;

G_D = 1/(Rs+Ld*s);
G_Q = 1/(Rs+Lq*s);

G_D = c2d(1/(Rs+Ld*s) * exp(-0.5*Ts*s),Ts,'zoh');
G_Q = c2d(1/(Rs+Lq*s) * exp(-0.5*Ts*s),Ts,'zoh');

%阶跃统计序列段
Start = 50;
Stop = 500;
%% D轴正向阶跃
data = load("D_Step_P.txt");
input = data(:,8);
output_rec = data(:,5);
t = Ts:Ts:Ts*(length(input));
output_sim = lsim(G_D,input(1:end),t);
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop),LineWidth=2);
title('D轴开环正向阶跃响应对比');
legend('实测输出','仿真输出');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('D轴开环正向阶跃相关系数：%.5f%%\n\n',r*100);

%% D轴负向阶跃
data = load("D_Step_N.txt");
input = data(:,8);
output_rec = data(:,5);
output_sim = lsim(G_D,input(1:end),t);
figure;
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop),LineWidth=2);
title('D轴开环负向阶跃响应对比');
legend('实测输出','仿真输出');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('D轴开环负向阶跃相关系数：%.5f%%\n\n',r*100);

%% Q轴正向阶跃
data = load("Q_Step_P.txt");
Uq = data(:,9);
theta_e = data(:,2);
Id = data(:,5);
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
output_rec = data(:,6);
output_sim = lsim(G_Q,input(1:end),t);
figure;
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop),LineWidth=2);
title('Q轴开环正向阶跃响应对比');
legend('实测输出','仿真输出');
grid on;
grid minor;
ylabel('电流/A');
xlabel('时间/s');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('Q轴开环正向阶跃相关系数：%.5f%%\n\n',r*100);

%% Q轴负向阶跃
data = load("Q_Step_N.txt");
Uq = data(:,9);
theta_e = data(:,2);
Id = data(:,5);
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
output_rec = data(:,6);
output_sim = lsim(G_Q,input(1:end),t);
figure;
plot(t(Start:Stop),output_rec(Start:Stop), ...
    t(Start:Stop),output_sim(Start:Stop),LineWidth=2);
grid on;
grid minor;
title('Q轴开环负向阶跃响应对比');
legend('实测输出','仿真输出');
ylabel('电流/A');
xlabel('时间/s');
cov_sim_rec = cov(output_rec(Start:Stop), output_sim(Start:Stop));
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(output_rec(Start:Stop)))* ...
    sqrt(var(output_sim(Start:Stop))));
fprintf('Q轴开环负向阶跃相关系数：%.5f%%\n\n',r*100);