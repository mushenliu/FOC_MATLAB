close all;
clear;
clc;
%% 数据处理函数
function [A,E,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,Y,n,k)
%输入参数定义
%X      输入变量矩阵
%Y      输出向量
%n      样本数
%k      待估计参数个数（自由度）
%输出参数定义
%A      估计参数向量
%E      残差向量
%Yhat   估计输出向量
%R2     决定系数
%R2_ADJ 调整决定系数
%Ld_hat Ld估计值
%Ld_sigma Ld标准差
%Lq_hat Lq估计值
%Lq_sigma Lq标准差

A = (X'*X) \ X' * Y;
Yhat = X * A;
E = Y - Yhat;
%残差平方和——无法解释的异变
SSE = sum(E.^2);
%总异变
SST = sum((Y - mean(Y)).^2);
%决定系数
R2 = 1 -SSE / SST;
%调整决定系数
R2_ADJ = 1 - (SSE / (n - k)) / (SST / (n - 1));
%残差分布估计E~N(0,σ^2)
sigma_hat_2_E = SSE / (n-k);
%参数向量的分布估计
sigma_hat_2_A = sigma_hat_2_E * inv(X'*X);
% Ld估计
Ld_hat = (A(1)-sqrt(A(2)^2+A(3)^2))/2;
Grad_Ld = [0.5,-0.5*A(2)/(sqrt(A(2)^2+A(3)^2)),-0.5*A(3)/(sqrt(A(2)^2+A(3)^2))];
Ld_sigma = sqrt(Grad_Ld * sigma_hat_2_A * Grad_Ld');
%Lq估计
Lq_hat = (A(1)+sqrt(A(2)^2+A(3)^2))/2;
Grad_Lq = [0.5,0.5*A(2)/(sqrt(A(2)^2+A(3)^2)),0.5*A(3)/(sqrt(A(2)^2+A(3)^2))];
Lq_sigma = sqrt(Grad_Lq * sigma_hat_2_A * Grad_Lq');

end
%% 数据提取
filename = '相间电感变化曲线.xlsx';
sheet = '原始数据';
data = readtable(filename, 'Sheet', sheet, 'Range', 'A3:G38', ...
    'VariableNamingRule', 'preserve');
% 指定变量名
data.Properties.VariableNames = {'theta_deg', 'AB_1k', 'BC_1k', 'AC_1k', ...
                                 'AB_10k', 'BC_10k', 'AC_10k'};
theta = data.theta_deg;
cos_2theta = cos(2*deg2rad(theta));
sin_2theta = sin(2*deg2rad(theta));
X = [ones(1,length(cos_2theta));cos_2theta';sin_2theta']';

L_AB_1k = data.AB_1k;
L_BC_1k = data.BC_1k;
L_AC_1k = data.AC_1k;
L_AB_10k = data.AB_10k;
L_BC_10k = data.BC_10k;
L_AC_10k = data.AC_10k;

%样本量
n = length(L_AC_10k);
%待估计参数个数（包含截距项）
k = 3;
%估测σ上下界
kp = 3;%上界
kn = 3;%下界

%结果记录
Ld_est = zeros(6,1);
Ld_sig = zeros(6,1);
Lq_est = zeros(6,1);
Lq_sig = zeros(6,1);
Y = zeros(6,length(theta));

%% 1kHz数据多元线性回归计算
fprintf('AB相@1kHz拟合结果\n');
[A,~,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,L_AB_1k,n,k);
fprintf(['拟合表达式：L=%.5f%+.5fcos2θ%+.5fsin2θ\n=' ...
    '%.5f%+.5fcos(2θ%+.5f)\n'], ...
    A(1),A(2),A(3),A(1), ...
    sqrt(A(2)^2+A(3)^2),rad2deg(atan2(-A(3),A(2))));
fprintf('决定系数R^2=%.5f\n',R2);
fprintf('调整决定系数R^2_ADJ=%.5f\n',R2_ADJ);
fprintf('Ld估计值：%.5fmH，标准差：%.5fmH\n',Ld_hat,Ld_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]\n',Ld_hat-3*Ld_sigma,Ld_hat+3*Ld_sigma);
fprintf('Lq估计值：%.5fmH，标准差：%.5fmH\n',Lq_hat,Lq_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]',Lq_hat-3*Lq_sigma,Lq_hat+3*Lq_sigma);
Y(1,:) = Yhat;
Ld_est(1) = Ld_hat;
Ld_sig(1) = Ld_sigma;
Lq_est(1) = Lq_hat;
Lq_sig(1) = Lq_sigma;

fprintf('\n\nBC相@1kHz拟合结果\n');
[A,~,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,L_BC_1k,n,k);
fprintf(['拟合表达式：L=%.5f%+.5fcos2θ%+.5fsin2θ\n=' ...
    '%.5f%+.5fcos(2θ%+.5f)\n'], ...
    A(1),A(2),A(3),A(1), ...
    sqrt(A(2)^2+A(3)^2),rad2deg(atan2(-A(3),A(2))));
fprintf('决定系数R^2=%.5f\n',R2);
fprintf('调整决定系数R^2_ADJ=%.5f\n',R2_ADJ);
fprintf('Ld估计值：%.5fmH，标准差：%.5fmH\n',Ld_hat,Ld_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]\n',Ld_hat-3*Ld_sigma,Ld_hat+3*Ld_sigma);
fprintf('Lq估计值：%.5fmH，标准差：%.5fmH\n',Lq_hat,Lq_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]',Lq_hat-3*Lq_sigma,Lq_hat+3*Lq_sigma);
Y(2,:) = Yhat;
Ld_est(2) = Ld_hat;
Ld_sig(2) = Ld_sigma;
Lq_est(2) = Lq_hat;
Lq_sig(2) = Lq_sigma;

fprintf('\n\nAC相@1kHz拟合结果\n');
[A,~,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,L_AC_1k,n,k);
fprintf(['拟合表达式：L=%.5f%+.5fcos2θ%+.5fsin2θ\n=' ...
    '%.5f%+.5fcos(2θ%+.5f)\n'], ...
    A(1),A(2),A(3),A(1), ...
    sqrt(A(2)^2+A(3)^2),rad2deg(atan2(-A(3),A(2))));
fprintf('决定系数R^2=%.5f\n',R2);
fprintf('调整决定系数R^2_ADJ=%.5f\n',R2_ADJ);
fprintf('Ld估计值：%.5fmH，标准差：%.5fmH\n',Ld_hat,Ld_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]\n',Ld_hat-3*Ld_sigma,Ld_hat+3*Ld_sigma);
fprintf('Lq估计值：%.5fmH，标准差：%.5fmH\n',Lq_hat,Lq_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]',Lq_hat-3*Lq_sigma,Lq_hat+3*Lq_sigma);
Y(3,:) = Yhat;
Ld_est(3) = Ld_hat;
Ld_sig(3) = Ld_sigma;
Lq_est(3) = Lq_hat;
Lq_sig(3) = Lq_sigma;

%% 10kHz数据多元线性回归计算
fprintf('\n\nAB相@10kHz拟合结果\n');
[A,~,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,L_AB_10k,n,k);
fprintf(['拟合表达式：L=%.5f%+.5fcos2θ%+.5fsin2θ\n=' ...
    '%.5f%+.5fcos(2θ%+.5f)\n'], ...
    A(1),A(2),A(3),A(1), ...
    sqrt(A(2)^2+A(3)^2),rad2deg(atan2(-A(3),A(2))));
fprintf('决定系数R^2=%.5f\n',R2);
fprintf('调整决定系数R^2_ADJ=%.5f\n',R2_ADJ);
fprintf('Ld估计值：%.5fmH，标准差：%.5fmH\n',Ld_hat,Ld_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]\n',Ld_hat-3*Ld_sigma,Ld_hat+3*Ld_sigma);
fprintf('Lq估计值：%.5fmH，标准差：%.5fmH\n',Lq_hat,Lq_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]',Lq_hat-3*Lq_sigma,Lq_hat+3*Lq_sigma);
Y(4,:) = Yhat;
Ld_est(4) = Ld_hat;
Ld_sig(4) = Ld_sigma;
Lq_est(4) = Lq_hat;
Lq_sig(4) = Lq_sigma;

fprintf('\n\nBC相@10kHz拟合结果\n');
[A,~,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,L_BC_10k,n,k);
fprintf(['拟合表达式：L=%.5f%+.5fcos2θ%+.5fsin2θ\n=' ...
    '%.5f%+.5fcos(2θ%+.5f)\n'], ...
    A(1),A(2),A(3),A(1), ...
    sqrt(A(2)^2+A(3)^2),rad2deg(atan2(-A(3),A(2))));
fprintf('决定系数R^2=%.5f\n',R2);
fprintf('调整决定系数R^2_ADJ=%.5f\n',R2_ADJ);
fprintf('Ld估计值：%.5fmH，标准差：%.5fmH\n',Ld_hat,Ld_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]\n',Ld_hat-3*Ld_sigma,Ld_hat+3*Ld_sigma);
fprintf('Lq估计值：%.5fmH，标准差：%.5fmH\n',Lq_hat,Lq_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]',Lq_hat-3*Lq_sigma,Lq_hat+3*Lq_sigma);
Y(5,:) = Yhat;
Ld_est(5) = Ld_hat;
Ld_sig(5) = Ld_sigma;
Lq_est(5) = Lq_hat;
Lq_sig(5) = Lq_sigma;

fprintf('\n\nAC相@10kHz拟合结果\n');
[A,~,Yhat,R2,R2_ADJ,Ld_hat,Ld_sigma,Lq_hat,Lq_sigma] = DataProcess (X,L_AC_10k,n,k);
fprintf(['拟合表达式：L=%.5f%+.5fcos2θ%+.5fsin2θ\n=' ...
    '%.5f%+.5fcos(2θ%+.5f)\n'], ...
    A(1),A(2),A(3),A(1), ...
    sqrt(A(2)^2+A(3)^2),rad2deg(atan2(-A(3),A(2))));
fprintf('决定系数R^2=%.5f\n',R2);
fprintf('调整决定系数R^2_ADJ=%.5f\n',R2_ADJ);
fprintf('Ld估计值：%.5fmH，标准差：%.5fmH\n',Ld_hat,Ld_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]\n',Ld_hat-3*Ld_sigma,Ld_hat+3*Ld_sigma);
fprintf('Lq估计值：%.5fmH，标准差：%.5fmH\n',Lq_hat,Lq_sigma);
fprintf('3σ的99.7%%置信度区间为：[%.5f,%.5f]',Lq_hat-3*Lq_sigma,Lq_hat+3*Lq_sigma);
Y(6,:) = Yhat;
Ld_est(6) = Ld_hat;
Ld_sig(6) = Ld_sigma;
Lq_est(6) = Lq_hat;
Lq_sig(6) = Lq_sigma;

%% 拟合曲线
figure(Name='1kHz线性回归曲线');
subplot(3,1,1);
plot(theta,L_AB_1k,theta,Y(1,:),LineWidth=2);
grid on;
grid minor;
ylabel('电感值/mH');
subtitle('AB')
legend('测量曲线','拟合曲线');
subplot(3,1,2);
plot(theta,L_BC_1k,theta,Y(2,:),LineWidth=2);
grid on;
grid minor;
ylabel('电感值/mH');
subtitle('BC')
legend('测量曲线','拟合曲线');
subplot(3,1,3);
plot(theta,L_AC_1k,theta,Y(3,:),LineWidth=2);
grid on;
grid minor;
ylabel('电感值/mH');
subtitle('AC')
legend('测量曲线','拟合曲线');
sgtitle('相间电感拟合结果（1kHz）');

figure(Name='10kHz线性回归曲线');
subplot(3,1,1);
plot(theta,L_AB_10k,theta,Y(4,:),LineWidth=2);
grid on;
grid minor;
ylabel('电感值/mH');
subtitle('AB')
legend('测量曲线','拟合曲线');
subplot(3,1,2);
plot(theta,L_BC_10k,theta,Y(5,:),LineWidth=2);
grid on;
grid minor;
ylabel('电感值/mH');
subtitle('BC')
legend('测量曲线','拟合曲线');
subplot(3,1,3);
plot(theta,L_AC_10k,theta,Y(6,:),LineWidth=2);
grid on;
grid minor;
ylabel('电感值/mH');
subtitle('AC')
legend('测量曲线','拟合曲线');
sgtitle('相间电感拟合结果（10kHz）');

%% 同步电感估测结果

% 分组索引
idx_1k = 1:3;   % AB, BC, AC @1kHz
idx_10k = 4:6;  % AB, BC, AC @10kHz

% 横坐标位置
x = 1:3;
labels = {'AB', 'BC', 'AC'};

figure(Name='置信区间');

subplot(2,2,1);
hold on;
% 误差棒（中心线、竖线、上下端线）
h = errorbar(x, Ld_est(idx_1k), kn * Ld_sig(idx_1k), kp * Ld_sig(idx_1k),...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410], ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'CapSize', 10);
% 额外加粗上下界水平线（使边界更醒目）
for i = 1:3
    % 上界
    plot([x(i)-0.15, x(i)+0.15], ...
         [Ld_est(idx_1k(i)) + kp*Ld_sig(idx_1k(i)), Ld_est(idx_1k(i)) + kp*Ld_sig(idx_1k(i))], ...
         'r-', 'LineWidth', 2.5);
    % 下界
    plot([x(i)-0.15, x(i)+0.15], ...
         [Ld_est(idx_1k(i)) - kn*Ld_sig(idx_1k(i)), Ld_est(idx_1k(i)) - kn*Ld_sig(idx_1k(i))], ...
         'r-', 'LineWidth', 2.5);
end
% 修饰
set(gca, 'XTick', x, 'XTickLabel', labels);
ylabel('Ld/mH');
title(sprintf('Ld@1kHz'));
grid on;
hold off;

% ---- 子图2: Ld @ 10kHz ----
subplot(2,2,2);
hold on;
h = errorbar(x, Ld_est(idx_10k),kn * Ld_sig(idx_10k),kp * Ld_sig(idx_10k), ...
    's', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410], ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'CapSize', 10);
for i = 1:3
    plot([x(i)-0.15, x(i)+0.15], ...
         [Ld_est(idx_10k(i)) + kp*Ld_sig(idx_10k(i)), Ld_est(idx_10k(i)) + kp*Ld_sig(idx_10k(i))], ...
         'r-', 'LineWidth', 2.5);
    plot([x(i)-0.15, x(i)+0.15], ...
         [Ld_est(idx_10k(i)) - kn*Ld_sig(idx_10k(i)), Ld_est(idx_10k(i)) - kn*Ld_sig(idx_10k(i))], ...
         'r-', 'LineWidth', 2.5);
end
set(gca, 'XTick', x, 'XTickLabel', labels);
ylabel('Ld/mH');
title(sprintf('Ld@10kHz'));
grid on;
hold off;

% ---- 子图3: Lq @ 1kHz ----
subplot(2,2,3);
hold on;
h = errorbar(x, Lq_est(idx_1k),kn * Lq_sig(idx_1k),kp * Lq_sig(idx_1k), ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'CapSize', 10);
for i = 1:3
    plot([x(i)-0.15, x(i)+0.15], ...
         [Lq_est(idx_1k(i)) + kp*Lq_sig(idx_1k(i)), Lq_est(idx_1k(i)) + kp*Lq_sig(idx_1k(i))], ...
         'b-', 'LineWidth', 2.5);
    plot([x(i)-0.15, x(i)+0.15], ...
         [Lq_est(idx_1k(i)) - kn*Lq_sig(idx_1k(i)), Lq_est(idx_1k(i)) - kn*Lq_sig(idx_1k(i))], ...
         'b-', 'LineWidth', 2.5);
end
set(gca, 'XTick', x, 'XTickLabel', labels);
ylabel('Lq/mH');
title(sprintf('Lq@1kHz'));
grid on;
hold off;

% ---- 子图4: Lq @ 10kHz ----
subplot(2,2,4);
hold on;
h = errorbar(x, Lq_est(idx_10k),kn * Lq_sig(idx_10k),kp * Lq_sig(idx_10k), ...
    's', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'CapSize', 10);
for i = 1:3
    plot([x(i)-0.15, x(i)+0.15], ...
         [Lq_est(idx_10k(i)) + kp*Lq_sig(idx_10k(i)), Lq_est(idx_10k(i)) + kp*Lq_sig(idx_10k(i))], ...
         'b-', 'LineWidth', 2.5);
    plot([x(i)-0.15, x(i)+0.15], ...
         [Lq_est(idx_10k(i)) - kn*Lq_sig(idx_10k(i)), Lq_est(idx_10k(i)) - kn*Lq_sig(idx_10k(i))], ...
         'b-', 'LineWidth', 2.5);
end
set(gca, 'XTick', x, 'XTickLabel', labels);
ylabel('Lq/mH');
title(sprintf('Lq@10kHz'));
grid on;
hold off;

sgtitle(sprintf('同步电感估计值及 -%dσ~%dσ 置信区间', kn,kp));

%% 对频率敏感性分析
fprintf(['\n\n10kHz三组Ld估计均值与1kHz三组Ld估计均值之差' ...
    '占整体六组Ld估计均值之比：%.5f%%'], ...
    (mean(Ld_est(1:3)) - mean(Ld_est(4:6))) / mean(Ld_est) * 100);
fprintf(['\n\n六组Ld估计极差' ...
    '占整体六组Ld估计均值之比：%.5f%%'], ...
    (max(Ld_est) - min(Ld_est)) / mean(Ld_est) * 100);
fprintf(['\n\n10kHz三组Lq估计均值与1kHz三组Lq估计均值之差' ...
    '占整体六组Ld估计均值之比：%.5f%%'], ...
    (mean(Lq_est(1:3)) - mean(Lq_est(4:6))) / mean(Lq_est) * 100);
fprintf(['\n\n六组Lq估计极差' ...
    '占整体六组Lq估计均值之比：%.5f%%'], ...
    (max(Lq_est) - min(Lq_est)) / mean(Lq_est) * 100);
fprintf(['\n\n六组Ld估计结果' ...
    '（均值±0.5极差）：%.5f mH~%.5f mH'], ...
    mean(Ld_est)-0.5*(max(Ld_est) - min(Ld_est)), ...
    mean(Ld_est)+0.5*(max(Ld_est) - min(Ld_est)));
fprintf(['\n\n六组Lq估计结果' ...
    '（均值±0.5极差）：%.5f mH~%.5f mH'], ...
    mean(Lq_est)-0.5*(max(Lq_est) - min(Lq_est)), ...
    mean(Lq_est)+0.5*(max(Lq_est) - min(Lq_est)));