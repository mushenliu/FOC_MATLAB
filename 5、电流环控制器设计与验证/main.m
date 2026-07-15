close all;
clear;
clc;
Ts = 1/(170000000 / 5665 / 2);
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.PhaseWrapping = 'on';
opts.Grid = 'on';
opts.XLim = [1, 1/(2*Ts)];

s=tf('s');
z=tf('z',Ts);
%阶跃统计序列段
D_Start = 990;
D_Stop = 1100;
Q_Start = 990;
Q_Stop = 1025;
%% 参数辨识结果
Rs = 5.70125;
Ld = 1.86152e-3;
Psi = 9.65548e-3;
Lq = 2.16828e-3;

G_D = c2d(1/(Rs+Ld*s) * exp(-1*Ts*s) * exp(-260e-9 * s) * 1/(120e-9 * s + 1), ...
    Ts,'zoh');
G_Q = c2d(1/(Rs+Lq*s) * exp(-1*Ts*s) * exp(-260e-9 * s) * 1/(120e-9 * s + 1), ...
    Ts,'zoh');
%% D轴频域矫正控制器设计
%相位补偿/°
phi_d = 60;
%补偿频率/Hz
f_d = 1000;
%计算频域矫正控制器
alpha_d = (1-sin(deg2rad(phi_d)))/(1+sin(deg2rad(phi_d)));
Td = 1/(2*pi*f_d*sqrt(alpha_d));
Kd = 1;
Cd = Kd * 1/s * (Td*s+1)/(alpha_d*Td*s+1);
Cd = c2d(Cd,Ts,'tustin');
[A,~,~]=bode(G_D*Cd,2*pi*f_d);
Kd = 1/A;
Cd = Cd * Kd;
figure(Name='D轴伯德图');
subplot(1,2,1);
bode(G_D,G_D*c2d(1/s,Ts,'tustin'),G_D*Cd,opts);
grid on;
grid minor;
subtitle('D轴电流环开环伯德图');
legend('被控对象','带积分环节的等价被控对象','带控制环节的开环伯德图');
Close_d = feedback(G_D*Cd,1);
subplot(1,2,2);
bode(Close_d,opts);
grid on;
grid minor;
subtitle('D轴电流环闭环伯德图');
figure(Name='D轴开闭环对比分析');
subplot(1,2,1);
step(G_D,Close_d);
grid on;
grid minor;
subtitle('D轴电流环单位阶跃响应');
legend('开环','闭环');
subplot(1,2,2);
pzmap(G_D,Close_d);
subtitle('零极点图');
legend('开环','闭环');
[Gm,Pm,Wcg,Wcp] = margin(G_D*Cd);
BW = bandwidth(Close_d);
fprintf('D轴电流环超前校正：\n')
fprintf('开环截止频率：%.5fHz\n',Wcp/(2*pi));
fprintf('相位裕度：%.5f°\n',Pm);
fprintf('相位穿越频率：%.5fHz\n',Wcg/(2*pi));
fprintf('幅值裕度：%.5fdB\n',20*log10(Gm));
fprintf('闭环带宽：%.5fHz\n',BW/(2*pi));

%% Q轴频域矫正控制器设计
%相位补偿/°
phi_q = 60;
%补偿频率/Hz
f_q = 1000;
%计算频域矫正控制器
alpha_q = (1-sin(deg2rad(phi_q)))/(1+sin(deg2rad(phi_q)));
Tq = 1/(2*pi*f_q*sqrt(alpha_q));
Kq = 1;
Cq = Kq * 1/s * (Tq*s+1)/(alpha_q*Tq*s+1);
Cq = c2d(Cq,Ts,'tustin');
[A,~,~]=bode(G_Q*Cq,2*pi*f_q);
Kq = 1/A;
Cq = Cq * Kq;
figure(Name='Q轴伯德图');
subplot(1,2,1);
bode(G_Q,G_Q*c2d(1/s,Ts,'tustin'),G_Q*Cq,opts);
grid on;
grid minor;
subtitle('Q轴电流环开环伯德图');
legend('被控对象','带积分环节的等价被控对象','带控制环节的开环伯德图');
Close_q = feedback(G_Q*Cq,1);
subplot(1,2,2);
bode(Close_q,opts);
grid on;
grid minor;
subtitle('Q轴电流环闭环伯德图');
figure(Name='Q轴开闭环对比分析');
subplot(1,2,1);
step(G_Q,Close_q);
grid on;
grid minor;
subtitle('Q轴电流环单位阶跃响应');
legend('开环','闭环');
subplot(1,2,2);
pzmap(G_Q,Close_q);
subtitle('零极点图');
legend('开环','闭环');
[Gm,Pm,Wcg,Wcp] = margin(G_Q*Cq);
BW = bandwidth(Close_q);
fprintf('Q轴电流环超前校正：\n')
fprintf('开环截止频率：%.5fHz\n',Wcp/(2*pi));
fprintf('相位裕度：%.5f°\n',Pm);
fprintf('相位穿越频率：%.5fHz\n',Wcg/(2*pi));
fprintf('幅值裕度：%.5fdB\n',20*log10(Gm));
fprintf('闭环带宽：%.5fHz\n',BW/(2*pi));

%% D轴闭环正向验证
data =  readmatrix("D_Close_P.txt", 'NumHeaderLines', 1);
input = data(end-2046:end,10);
output_rec = data(end-2046:end,5);
Iq = data(end-2046:end,6);
theta_e = data(end-2046:end,2);
we = zeros(length(theta_e),1);
we(1) = 0;
for i=2:length(we)
    we(i)=(theta_e(i)-theta_e(i-1))/Ts * 2 * pi / 360;
    if(abs(theta_e(i)-theta_e(i-1))>50)
        we(i)=we(i-1);
    end
end
t = Ts:Ts:Ts*(length(input));
output_sim = lsim(Close_d,input,t);
figure(Name='D轴正向');
plot(t(D_Start:D_Stop),output_rec(D_Start:D_Stop), ...
    t(D_Start:D_Stop),output_sim(D_Start:D_Stop), ...
    t(D_Start:D_Stop),input(D_Start:D_Stop),LineWidth=2);
title('D轴正向方波给定的闭环响应对比');
legend('实测输出','仿真输出','给定输入');
grid on;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nD轴正向方波响应分析：\n')
Signal_Analyse(output_rec(D_Start:D_Stop),output_sim(D_Start:D_Stop));

%% D轴闭环负向验证
data =  readmatrix("D_Close_N.txt", 'NumHeaderLines', 1);
input = data(end-2046:end,10);
output_rec = data(end-2046:end,5);
t = Ts:Ts:Ts*(length(input));
output_sim = lsim(Close_d,input,t);
figure(Name='D轴负向');
plot(t(D_Start:D_Stop),output_rec(D_Start:D_Stop), ...
    t(D_Start:D_Stop),output_sim(D_Start:D_Stop), ...
    t(D_Start:D_Stop),input(D_Start:D_Stop),LineWidth=2);
title('D轴负向方波给定的闭环响应对比');
legend('实测输出','仿真输出','给定输入');
grid on;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nD轴负向方波响应分析：\n')
Signal_Analyse(output_rec(D_Start:D_Stop),output_sim(D_Start:D_Stop));

%% Q轴闭环正向验证
data =  readmatrix("Q_Close_P.txt", 'NumHeaderLines', 1);
input = data(end-2046:end,11);
output_rec = data(end-2046:end,6);
t = Ts:Ts:Ts*(length(input));
output_sim = lsim(Close_q,input,t);
figure(Name='Q轴正向');
plot(t(Q_Start:Q_Stop),output_rec(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),input(Q_Start:Q_Stop),LineWidth=2);
title('Q轴正向方波给定的闭环响应对比');
legend('实测输出','仿真输出','给定输入');
grid on;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nQ轴正向方波响应分析：\n')
Signal_Analyse(output_rec(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop));

%% Q轴闭环负向验证
data =  readmatrix("Q_Close_N.txt", 'NumHeaderLines', 1);
input = data(end-2046:end,11);
output_rec = data(end-2046:end,6);
t = Ts:Ts:Ts*(length(input));
output_sim = lsim(Close_q,input,t);
figure(Name='Q轴负向');
plot(t(Q_Start:Q_Stop),output_rec(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop), ...
    t(Q_Start:Q_Stop),input(Q_Start:Q_Stop),LineWidth=2);
title('Q轴负向方波给定的闭环响应对比');
legend('实测输出','仿真输出','给定输入');
grid on;
ylabel('电流/A');
xlabel('时间/s');
fprintf('\nQ轴负向方波响应分析：\n')
Signal_Analyse(output_rec(Q_Start:Q_Stop),output_sim(Q_Start:Q_Stop));