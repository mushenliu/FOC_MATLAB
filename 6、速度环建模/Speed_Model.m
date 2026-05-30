close all;
clear;
clc;  
%% HK辨识函数定义
function [A, B, C, D, sys_ss, sys_tf] = HK (input, output, k, m, Ts)
%输入参数定义
%input    输入信号（已经去除暂态）
%output   输出信号（已经去除暂态）
%k        汉克尔矩阵阶次
%m        系统阶次
%Ts       采样时间 
%输出参数定义
%A        状态矩阵
%B        输入矩阵
%C        输出矩阵
%D        馈通矩阵
%sys_ss   状态空间
%sys_tf   传递函数

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';

%一、系统脉冲响应序列求解
input = input - mean(input);
output = output - mean(output);
%去除直流分量
n = length(input);
Ruu = xcorr(input, input)';
Ryu = xcorr(output, input)';
figure;
hold on;
title('输入输出互相关函数');
i=-(n-1):1:n-1;
plot(i,Ryu);
figure;
hold on;
title('输入自相关函数函数');
plot(i,Ruu);
%xcorr计算得到的是列向量，需要转置得到行向量，且原本的[-(n-1),n-1]索引范围会变成[1,2n-1]
%即全部加n，并且互相关计算时要把输出信号前置，输入信号后置（xcorr存在顺序问题）
for i = 1:n%构造自相关托普利茨矩阵
    ruu(:,i) = Ruu(n+1-i:2*n-i);
end
ryu = Ryu(n:2*n-1);%构造互相关列向量
g = inv(ruu)*ryu';%系统脉冲响应序列，序号从1开始，因此g(k)其实在k+1的位置上
figure;
plot(g);
title('脉冲响应序列');

%二、系统阶次分析
for i = 1:k%构造汉克尔矩阵
    H(:,i) = g(i+1:i+k);
end
[U, S, V]=svd(H);%奇异值分解
sigma = S.^(0.5);
S = 1/S(1,1) .* S;
figure;
title('汉克尔矩阵奇异值');
hold on;
for i=1:k
    stem(i,S(i,i));
end

%三、系统模型求解
for i = 1:k%构造移位汉克尔矩阵
    H1(:,i) = g(i+2:i+k+1);
end
%输入矩阵求解
BO = sigma * V';
B = BO(:,1);
B = B(1:m,:);
%输出矩阵求解
CO = U * sigma;
C = CO(1,:);
C = C(:,1:m);
%状态矩阵求解
AO = inv(CO)*H1*inv(BO);
A = AO(1:m,1:m);
%馈通矩阵求解
D = g(1);
sys_ss = ss(A,B,C,D,Ts);
sys_tf = tf(sys_ss);
end
%% 数据处理函数定义
function [sys_tf,Gest] = DataProcess (input, output, k, m, Ts)
%输入参数定义
%input    输入信号（已经去除暂态）
%output   输出信号（已经去除暂态）
%k        汉克尔矩阵阶次
%m        系统阶次
%Ts       采样时间 
%输出参数定义
%sys_tf   传递函数
%Gest     频率响应对象

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';    
opts.Grid = 'on';

n=length(input);
%频点向量
f = 1 / Ts ;
F = 1:1:(n); 
F =  F * f / (n);
w = 2*pi*F;

[~, ~, ~, ~, ~, sys_tf] = HK (input, output, k, m, Ts);

%估计系统的频率特性
%双重循环求解离散非周期信号的傅里叶变化
Ruu = xcorr(input,input)';
Ryu = xcorr(output,input)';
Suu = zeros(1,length(w));
for k = 1:n
    for j = 1:n
        Suu(k) = Suu(k) + Ruu(j+n-1) * exp(-1i * 2 * pi * k * j / (n));
    end
end
Syu = zeros(1,length(w));
for k = 1:n
    for j = 1:n
        Syu(k) = Syu(k) + Ryu(j+n-1) * exp(-1i * 2 * pi * k * j / (n));
    end
end
H_est = Syu./Suu;
Gest = frd(H_est,w);

figure;
bode(sys_tf,Gest,opts);
legend("解析解",'数值解');

[Ags,Pgs] = bode(Gest,w);
Ags = 20*log10(squeeze(Ags));
Pgs = squeeze(Pgs);
Pgs = mod(Pgs + 180, 360) - 180;

[Ahk,Phk] = bode(sys_tf,w);
Ahk = 20*log10(squeeze(Ahk));
Phk = squeeze(Phk);
Phk = mod(Phk + 180, 360) - 180;

cov_input_output = cov(Ahk, Ags);
cov_value=cov_input_output(1,2);
r=cov_value/(sqrt(var(Ahk))*sqrt(var(Ags)));
fprintf('幅频相关系数：%.5f%%\n',r*100);
cov_input_output = cov(Phk, Pgs);
cov_value=cov_input_output(1,2);
r=cov_value/(sqrt(var(Phk))*sqrt(var(Pgs)));
fprintf('相频相关系数：%.5f%%\n',r*100);

end
%% 速度环辨识
%Hankel矩阵阶次
k = 10;
%系统阶次
m = 2;
%采样时间
Ts = 1/1000;
s = tf('s');
z=tf('z',Ts);

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';     

data = load("Speed.txt");
input = data(:,11);
output = data(:,3);
[sys_tf,Gest] = DataProcess(input,output,k,m,Ts);
sprintf('传递函数：');
sys_tf
figure;
pzmap(sys_tf);
%% 速度环控制器设计
% %相位补偿/°
% phi = 50;
% %补偿频率/Hz
% f = 0.1;
% %补偿分贝数
% A = 28;
% %计算频域矫正控制器
% alpha_d = (1-sin(deg2rad(phi)))/(1+sin(deg2rad(phi)));
% Td = 1/(2*pi*f*sqrt(alpha_d));
% Kd = db2mag(A);
% Cd =Kd * (Td*s+1)/(alpha_d*Td*s+1);
% Cd = c2d(Cd,Ts,'tustin') * (z^2 - 1.672 *z + 0.7639)/((z-0.0382)*(z-0.5));
% figure;
% bode(sys_tf,sys_tf*Cd,opts);
% grid on;
% title('速度环开环伯德图');
% legend('被控对象','频域矫正');
% Close_s = feedback(sys_tf*Cd,1);
% figure;
% bode(Close_s,opts);
% grid on;
% title('速度环闭环伯德图');
% figure;
% step(sys_tf,Close_s);
% grid on;
% title('速度环单位阶跃响应');
% legend('被控对象','频域矫正');
% Cd

% C_PI = 0.001 + 0.005/s;
% C_PI = c2d(C_PI,Ts,'tustin');
% figure;
% bode(sys_tf,sys_tf*C_PI,opts);
% grid on;
% title('速度环开环伯德图');
% legend('被控对象','频域矫正');
% Close_s = feedback(sys_tf*C_PI,1);
% figure;
% bode(Close_s,opts);
% grid on;
% title('速度环闭环伯德图');
% figure;
% step(sys_tf,Close_s);
% grid on;
% title('速度环单位阶跃响应');
% legend('被控对象','频域矫正');
% C_PI