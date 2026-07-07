% Ho-Kalman辨识和辨识结果量化分析函数
function [sys_tf,Gest,r_AF,r_PF,BFT] = DataProcess (input, output, k, m, Ts)
%输入参数定义
%input    输入信号（已经去除暂态）
%output   输出信号（已经去除暂态）
%k        汉克尔矩阵阶次
%m        系统阶次
%Ts       采样时间
%输出参数定义
%sys_tf   传递函数
%Gest     频率响应对象
%r_AF     幅频相关系数
%r_PF     相频相关系数
%BFT      拟合优度

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
Ruu = xcorr(input,input)';
Ryu = xcorr(output,input)';
Suu = fft(Ruu(n:end));
Syu = fft(Ryu(n:end));
H_est = Syu./Suu;
Gest = frd(H_est,w);

figure(Name='频域对比');
bode(sys_tf,Gest,opts);
legend("解析解",'数值解');
grid on;
grid minor;

[Ags,Pgs] = bode(Gest,w);
Hgs = Ags .* exp(1j * deg2rad(Pgs));
Hgs = squeeze(Hgs);
Ags = squeeze(Ags);
Pgs = squeeze(Pgs);

[Ahk,Phk] = bode(sys_tf,w);
Hhk = Ahk .* exp(1j * deg2rad(Phk));
Hhk = squeeze(Hhk);
Ahk = squeeze(Ahk);
Phk = squeeze(Phk);

cov_input_output = cov(Ahk, Ags);
cov_value=cov_input_output(1,2);
r_AF=cov_value/(sqrt(var(Ahk))*sqrt(var(Ags)));
fprintf('幅频相关系数：%.5f%%\n',r_AF*100);
cov_input_output = cov(Phk, Pgs);
cov_value=cov_input_output(1,2);
r_PF=cov_value/(sqrt(var(Phk))*sqrt(var(Pgs)));
fprintf('相频相关系数：%.5f%%\n',r_PF*100);
BFT = (1 - norm(Hgs - Hhk) / norm(Hgs - mean(Hgs)));
fprintf('复频率响应拟合优度%.5f%%\n',BFT*100);

end
