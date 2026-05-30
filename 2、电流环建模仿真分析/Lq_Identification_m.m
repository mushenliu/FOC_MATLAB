close all;
clear;
clc;
%% 处理函数定义
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

%一、系统脉冲响应序列求解
input = input - mean(input);
output = output - mean(output);
% %去除直流分量
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

function [sys_tf,Gest,r_AF,r_PF] = DataProcess (input, output, k, m, Ts)
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
r_AF=cov_value/(sqrt(var(Ahk))*sqrt(var(Ags)))*100;
fprintf('幅频相关系数：%.5f%%\n',r_AF);
cov_input_output = cov(Phk, Pgs);
cov_value=cov_input_output(1,2);
r_PF=cov_value/(sqrt(var(Phk))*sqrt(var(Pgs)))*100;
fprintf('相频相关系数：%.5f%%\n',r_PF);

end
%% 仿真参数定义
Ld = 2e-3;
Lq = 2.5e-3;
Rs = 8;
Psi = 0.01;
np = 14;
Udc = 12;
J = 0.00062;
f = 15e3;
Ts = 1/f;
%伪随机信号级数
N = 12;
%汉克尔矩阵阶次
k = 20;
%系统阶次
m = 2;
%伪随机幅值
A = 1;
%PWM延迟拍数
Beat=1;

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';
%阶跃持续时间
Step = 2500;

%% 单工作点辨识
DC = 0;
fprintf('工作点%.0fV：\n',DC);
%生成输入信号
input_step = ones(Step,1)*DC;
t_step = (0:Ts:(Step-1)*Ts)';

input_prbs = idinput(2^N-1,'prbs');
for i=1:length(input_prbs)
    if(input_prbs(i)==1)
        input_prbs(i)=A+DC;
    else
        input_prbs(i)=-A+DC;
    end
end
t_prbs = (Step*Ts : Ts : (Step + 2^N - 2)*Ts)';

t_total = [t_step; t_prbs];
input_total = [input_step; input_prbs];
Uq = [t_total, input_total];

out = sim("Lq_Identification_sim.slx");
Id = out.Id.Data(Step + 1:end);
we = out.we.Data(Step + 1:end);
Uq_ac = input_prbs - we .* (Ld * Id + Psi);
output = out.Iq.Data(Step + 1:end);
[sys_tf_0,Gest_0] = DataProcess(Uq_ac,output,k,m,Ts);
Lq_0 = -Rs*Ts./log((pole(sys_tf_0)));
fprintf('辨识结果：%.5f mH\n相对偏差：%.5f%%\n\n', ...
    Lq_0(1)*1000, (Lq_0(1)-Lq)/Lq*100);

%% 多工作点辨识
% Lq_data = zeros(length(-5:0.1:5),1);
% Error = zeros(length(-5:0.1:5),1);
% AF = zeros(length(-5:0.1:5),1);
% PF = zeros(length(-5:0.1:5),1);
% point=1;
% U=-5:0.1:5;
% for i=-5:0.1:5
%     DC = i;
%     fprintf('工作点%.1fV：\n',DC);
%     %生成输入信号
%     input_step = ones(Step,1)*DC;
%     t_step = (0:Ts:(Step-1)*Ts)';
%     input_prbs = idinput(2^N-1,'prbs');
%         for j=1:length(input_prbs)
%             if(input_prbs(j)==1)
%                 input_prbs(j)=A+DC;
%             else
%                 input_prbs(j)=-A+DC;
%             end
%         end
%     t_prbs = (Step*Ts : Ts : (Step + 2^N - 2)*Ts)';
%     t_total = [t_step; t_prbs];
%     input_total = [input_step; input_prbs];
%     Uq = [t_total, input_total];
%     out = sim("Lq_Identification_sim.slx");
%     Id = out.Id.Data(Step + 1:end);
%     we = out.we.Data(Step + 1:end);
%     Uq_ac = input_prbs - we .* (Ld * Id + Psi);
%     output = out.Iq.Data(Step + 1:end);
%     [sys_tf,Gest,r_AF,r_PF] = DataProcess(Uq_ac,output,k,m,Ts);
%     Lq_o = -Rs*Ts./log((pole(sys_tf)));
%     fprintf('辨识结果：%.5f mH\n相对偏差：%.5f%%\n\n', ...
%         Lq_o(1)*1000, (Lq_o(1)-Lq)/Lq*100);
%     Lq_data(point,1)=Lq_o(1);
%     Error(point,1)=(Lq_o(1)-Lq)/Lq*100;
%     AF(point,1)=r_AF;
%     PF(point,1)=r_PF;
%     point = point+1;
%     close all;
% end
% figure;
% plot(U,AF,U,PF);
% title('幅频-相频相关系数');
% xlabel('D轴工作点电压/V');
% ylabel('相关系数/%');
% legend('幅频相关系数','相频相关系数');
% figure;
% plot(U,Error);
% title('不同工作点电感辨识相对误差');
% xlabel('Q轴工作点电压/V');
% ylabel('相对误差/%');