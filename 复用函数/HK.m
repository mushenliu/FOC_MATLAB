% Ho-Kalman辨识函数
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

figure(Name='Ho-Kalman辨识分析')
sgtitle('Ho-Kalman辨识分析');

%一、系统脉冲响应序列求解
input = input - mean(input);
output = output - mean(output);
%去除直流分量
n = length(input);
Ruu = xcorr(input, input)';
Ryu = xcorr(output, input)';
subplot(2,2,1);
hold on;
subtitle('输入输出互相关函数');
i=-(n-1):1:n-1;
plot(i,Ryu);
grid on;
grid minor;
subplot(2,2,2);
hold on;
subtitle('输入自相关函数函数');
plot(i,Ruu);
grid on;
grid minor;
%xcorr计算得到的是列向量，需要转置得到行向量，且原本的[-(n-1),n-1]索引范围会变成[1,2n-1]
%即全部加n，并且互相关计算时要把输出信号前置，输入信号后置（xcorr存在顺序问题）
for i = 1:n%构造自相关托普利茨矩阵
    ruu(:,i) = Ruu(n+1-i:2*n-i);
end
ryu = Ryu(n:2*n-1);%构造互相关列向量
g = ruu \ ryu';%系统脉冲响应序列，序号从1开始，因此g(k)其实在k+1的位置上
subplot(2,2,3);
plot(g);
subtitle('脉冲响应序列');
grid on;
grid minor;

%二、系统阶次分析
for i = 1:k%构造汉克尔矩阵
    H(:,i) = g(i+1:i+k);
end
[U, S, V]=svd(H);%奇异值分解
sigma = S.^(0.5);
S = 1/S(1,1) .* S;
subplot(2,2,4);
subtitle('汉克尔矩阵奇异值');
hold on;
for i=1:k
    stem(i,S(i,i));
end
grid on;
grid minor;

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
AO = CO \ H1 / BO;
A = AO(1:m,1:m);
%馈通矩阵求解
D = g(1);
sys_ss = ss(A,B,C,D,Ts);
sys_tf = tf(sys_ss);
end
