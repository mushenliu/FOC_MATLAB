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
% figure;
% hold on;
% title('输入输出互相关函数');
% i=-(n-1):1:n-1;
% plot(i,Ryu);
% figure;
% hold on;
% title('输入自相关函数函数');
% plot(i,Ruu);
%xcorr计算得到的是列向量，需要转置得到行向量，且原本的[-(n-1),n-1]索引范围会变成[1,2n-1]
%即全部加n，并且互相关计算时要把输出信号前置，输入信号后置（xcorr存在顺序问题）
for i = 1:n%构造自相关托普利茨矩阵
    ruu(:,i) = Ruu(n+1-i:2*n-i);
end
ryu = Ryu(n:2*n-1);%构造互相关列向量
g = inv(ruu)*ryu';%系统脉冲响应序列，序号从1开始，因此g(k)其实在k+1的位置上
% figure;
% plot(g);
% title('脉冲响应序列');

%二、系统阶次分析
for i = 1:k%构造汉克尔矩阵
    H(:,i) = g(i+1:i+k);
end
[U, S, V]=svd(H);%奇异值分解
sigma = S.^(0.5);
S = 1/S(1,1) .* S;
% figure;
% title('汉克尔矩阵奇异值');
% hold on;
% for i=1:k
%     stem(i,S(i,i));
% end

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
opts.XLim = [1, 1/(2*Ts)];

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
grid minor;

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

data = load("Ld_0.txt");
fprintf('工作点0V：\n');
input_p0 = data(:,8);
output_p0 = data(:,5) - mean(data(:,5));
[sys_tf_p0,Gest_p0] = DataProcess(input_p0,output_p0,k,m,Ts);
Ld_0 = -Rs_O*Ts./log((pole(sys_tf_p0)));
fprintf('辨识结果：%.5f mH\n\n',Ld_0(1)*1000);

data = load("Ld_2.txt");
fprintf('工作点2V：\n');
input_p2 = data(:,8) - 2;
output_p2 = data(:,5) - mean(data(:,5));
[sys_tf_p2,Gest_p2] = DataProcess(input_p2,output_p2,k,m,Ts);
Ld_p2 = -Rs_O*Ts./log((pole(sys_tf_p2)));
fprintf('辨识结果：%.5f mH\n\n',Ld_p2(1)*1000);

data = load("Ld_4.txt");
fprintf('工作点4V：\n');
input_p4 = data(:,8) - 4;
output_p4 = data(:,5) - mean(data(:,5));
[sys_tf_p4,Gest_p4] = DataProcess(input_p4,output_p4,k,m,Ts);
Ld_p4 = -Rs_O*Ts./log((pole(sys_tf_p4)));
fprintf('辨识结果：%.5f mH\n\n',Ld_p4(1)*1000);

data = load("Ld_-2.txt");
fprintf('工作点-2V：\n');
input_n2 = data(:,8) + 2;
output_n2 = data(:,5) - mean(data(:,5));
[sys_tf_n2,Gest_n2] = DataProcess(input_n2,output_n2,k,m,Ts);
Ld_n2 = -Rs_O*Ts./log((pole(sys_tf_n2)));
fprintf('辨识结果：%.5f mH\n\n',Ld_n2(1)*1000);

data = load("Ld_-4.txt");
fprintf('工作点-4V：\n');
input_n4 = data(:,8) + 4;
output_n4 = data(:,5) - mean(data(:,5));
[sys_tf_n4,Gest_n4] = DataProcess(input_n4,output_n4,k,m,Ts);
Ld_n4 = -Rs_O*Ts./log((pole(sys_tf_n4)));
fprintf('辨识结果：%.5f mH\n\n',Ld_n4(1)*1000);

figure;
pzmap(sys_tf_p0,sys_tf_p2,sys_tf_p4,sys_tf_n2,sys_tf_n4);
title('D轴传递函数零极点图');
grid on;
grid minor;
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');
figure;
bode(sys_tf_p0,sys_tf_p2,sys_tf_p4,sys_tf_n2,sys_tf_n4,opts);
title('D轴传递函数频率特性曲线');
grid on;
grid minor;
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');