close all;
clear;
clc;
%% HK参数设置
%Hankel矩阵阶次
k = 15;
%系统阶次
m = 2;
%采样时间
Ts = 1/(170000000 / 5665 / 2);
%滤波器系数
fi = 0.1;

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';
opts.XLim = [1, 1/(2*Ts)];

Rs_O = 5.70125;
Ld_O = 1.86152e-3;
Psi_O = 9.65548e-3;
Lq_O = 2.16828e-3;

%% 测量
data = readmatrix("Lq_0.txt", 'NumHeaderLines', 1);
fprintf('工作点0V：\n');
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
input_0 = Uq - we_lp .* (Ld_O * Id + Psi_O);
output_0 = data(end-2046:end,6);
[sys_tf_0,Gest_0] = DataProcess(input_0,output_0,k,m,Ts);
Lq_0 = -Rs_O*Ts./log((pole(sys_tf_0)));
fprintf('辨识结果：%.5f mH\n\n',Lq_0(1)*1000);

data = readmatrix("Lq_2.txt", 'NumHeaderLines', 1);
fprintf('工作点2V：\n');
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
input_p2 = Uq - we_lp .* (Ld_O * Id + Psi_O);
output_p2 = data(end-2046:end,6);
[sys_tf_p2,Gest_p2] = DataProcess(input_p2,output_p2,k,m,Ts);
Lq_p2 = -Rs_O*Ts./log((pole(sys_tf_p2)));
fprintf('辨识结果：%.5f mH\n\n',Lq_p2(1)*1000);

data = readmatrix("Lq_4.txt", 'NumHeaderLines', 1);
fprintf('工作点4V：\n');
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
input_p4 = Uq - we_lp .* (Ld_O * Id + Psi_O);
output_p4 = data(end-2046:end,6);
[sys_tf_p4,Gest_p4] = DataProcess(input_p4,output_p4,k,m,Ts);
Lq_p4 = -Rs_O*Ts./log((pole(sys_tf_p4)));
fprintf('辨识结果：%.5f mH\n\n',Lq_p4(1)*1000);

data = readmatrix("Lq_-2.txt", 'NumHeaderLines', 1);
fprintf('工作点-2V：\n');
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
input_n2 = Uq - we_lp .* (Ld_O * Id + Psi_O);
output_n2 = data(end-2046:end,6);
[sys_tf_n2,Gest_n2] = DataProcess(input_n2,output_n2,k,m,Ts);
Lq_n2 = -Rs_O*Ts./log((pole(sys_tf_n2)));
fprintf('辨识结果：%.5f mH\n\n',Lq_n2(1)*1000);

data = readmatrix("Lq_-4.txt", 'NumHeaderLines', 1);
fprintf('工作点-4V：\n');
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
input_n4 = Uq - we_lp .* (Ld_O * Id + Psi_O);
output_n4 = data(end-2046:end,6);
[sys_tf_n4,Gest_n4] = DataProcess(input_n4,output_n4,k,m,Ts);
Lq_n4 = -Rs_O*Ts./log((pole(sys_tf_n4)));
fprintf('辨识结果：%.5f mH\n\n',Lq_n4(1)*1000);

figure(Name='零极点图');
pzmap(sys_tf_0,sys_tf_p2,sys_tf_p4,sys_tf_n2,sys_tf_n4);
title('Q轴传递函数零极点图');
grid on;
grid minor;
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');
figure(Name='频率特性曲线');
bode(sys_tf_0,sys_tf_p2,sys_tf_p4,sys_tf_n2,sys_tf_n4,opts);
title('Q轴传递函数频率特性曲线');
grid on;
grid minor;
legend('工作点0','工作点2','工作点4','工作点-2','工作点-4');

%% 临时测试数据
% data = readmatrix("Lq_test.txt", 'NumHeaderLines', 1);
% Uq = data(end-2046:end,9);
% fprintf('工作点%.1f V：\n',mean(Uq));
% theta_e = data(end-2046:end,2);
% Id = data(end-2046:end,5);
% we = zeros(length(Id),1);
% we(1) = 0;
% for i=2:length(we)
%     we(i)=(theta_e(i)-theta_e(i-1))/Ts * 2 * pi / 360;
%     if(abs(theta_e(i)-theta_e(i-1))>50)
%         we(i)=we(i-1);
%     end
% end
% we_lp = lowpass(we,fi);
% input = Uq - we_lp .* (Ld_O * Id + Psi_O);
% output = data(end-2046:end,6);
% [sys_tf,Gest] = DataProcess(Uq,output,k,m,Ts);
% Lq = -Rs_O*Ts./log((pole(sys_tf)));
% fprintf('辨识结果：%.5f mH\n\n',Lq(1)*1000);
