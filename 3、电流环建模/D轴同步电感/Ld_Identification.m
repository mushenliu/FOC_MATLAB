close all;
clear;
clc;
%% HK参数设置
%Hankel矩阵阶次
k = 8;
%系统阶次
m = 2;
%采样时间
Ts = 1/(170000000 / 5665 / 2);

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';     

Rs_O = 5.70125;
%Ld = -Rs_O*Ts/log(pole)

%% 数据

data = readmatrix("Ld_0.txt", 'NumHeaderLines', 1);
fprintf('工作点0V：\n');
input_0 = data(end-2046:end,8);
output_0 = data(end-2046:end,5) - mean(data(end-2046:end,5));
[sys_tf_0,Gest_0] = DataProcess(input_0,output_0,k,m,Ts);
Ld_0 = -Rs_O*Ts./log((pole(sys_tf_0)));
fprintf('辨识结果：%.5f mH\n\n',Ld_0(1)*1000);

data = readmatrix("Ld_2.txt", 'NumHeaderLines', 1);
fprintf('工作点2V：\n');
input_p2 = data(end-2046:end,8) - 2;
output_p2 = data(end-2046:end,5) - mean(data(end-2046:end,5));
[sys_tf_p2,Gest_p2] = DataProcess(input_p2,output_p2,k,m,Ts);
Ld_p2 = -Rs_O*Ts./log((pole(sys_tf_p2)));
fprintf('辨识结果：%.5f mH\n\n',Ld_p2(1)*1000);

data = readmatrix("Ld_4.txt", 'NumHeaderLines', 1);
fprintf('工作点4V：\n');
input_p4 = data(end-2046:end,8) - 4;
output_p4 = data(end-2046:end,5) - mean(data(end-2046:end,5));
[sys_tf_p4,Gest_p4] = DataProcess(input_p4,output_p4,k,m,Ts);
Ld_p4 = -Rs_O*Ts./log((pole(sys_tf_p4)));
fprintf('辨识结果：%.5f mH\n\n',Ld_p4(1)*1000);

data = readmatrix("Ld_-2.txt", 'NumHeaderLines', 1);
fprintf('工作点-2V：\n');
input_n2 = data(end-2046:end,8) + 2;
output_n2 = data(end-2046:end,5) - mean(data(end-2046:end,5));
[sys_tf_n2,Gest_n2] = DataProcess(input_n2,output_n2,k,m,Ts);
Ld_n2 = -Rs_O*Ts./log((pole(sys_tf_n2)));
fprintf('辨识结果：%.5f mH\n\n',Ld_n2(1)*1000);

data = readmatrix("Ld_-4.txt", 'NumHeaderLines', 1);
fprintf('工作点-4V：\n');
input_n4 = data(end-2046:end,8) + 4;
output_n4 = data(end-2046:end,5) - mean(data(end-2046:end,5));
[sys_tf_n4,Gest_n4] = DataProcess(input_n4,output_n4,k,m,Ts);
Ld_n4 = -Rs_O*Ts./log((pole(sys_tf_n4)));
fprintf('辨识结果：%.5f mH\n\n',Ld_n4(1)*1000);

figure(Name='零极点图');
pzmap(sys_tf_0,sys_tf_p2,sys_tf_p4,sys_tf_n2,sys_tf_n4);
title('D轴传递函数零极点图');
grid on;
grid minor;
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');
figure(Name='频率特性曲线');
bode(sys_tf_0,sys_tf_p2,sys_tf_p4,sys_tf_n2,sys_tf_n4,opts);
title('D轴传递函数频率特性曲线');
grid on;
grid minor;
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');
figure(Name='相干函数对比');
[gamma_p0, F] = mscohere(input_0,output_0, hamming(256), 32, 1024, 1/Ts);
gamma_p0 = 20*log10(gamma_p0);
[gamma_p2, ~] = mscohere(input_p2,output_p2, hamming(256), 32, 1024, 1/Ts);
gamma_p2 = 20*log10(gamma_p2);
[gamma_p4, ~] = mscohere(input_p4,output_p4, hamming(256), 32, 1024, 1/Ts);
gamma_p4 = 20*log10(gamma_p4);
[gamma_n2, ~] = mscohere(input_n2,output_n2, hamming(256), 32, 1024, 1/Ts);
gamma_n2 = 20*log10(gamma_n2);
[gamma_n4, ~] = mscohere(input_n4,output_n4, hamming(256), 32, 1024, 1/Ts);
gamma_n4 = 20*log10(gamma_n4);
semilogx(F,gamma_p0,F,gamma_p2,F,gamma_p4,F,gamma_n2,F,gamma_n4);
grid on;
grid minor;
xlabel('频率/Hz');
ylabel('分贝/dB');
ylim([-0.5,0]);
xlim([100,3000]);
subtitle('输入输出相干函数');
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');

%% 临时测试数据
% data = readmatrix("Ld_test.txt", 'NumHeaderLines', 1);
% input = data(end-2046:end,8);
% fprintf('工作点%.1f V：\n',mean(input));
% input = data(end-2046:end,8) - mean(input);
% output = data(end-2046:end,5) - mean(data(end-2046:end,5));
% [sys_tf,Gest] = DataProcess(input,output,k,m,Ts);
% Ld = -Rs_O*Ts./log((pole(sys_tf)));
% fprintf('辨识结果：%.5f mH\n\n',Ld(1)*1000);