close all;
clear;
clc;
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

%理想与非理想仿真切换（≥0为非理想，＜0为理想）
Switch = -1;

%伪随机信号级数
N = 11;
%汉克尔矩阵阶次
k = 8;
%系统阶次
if Switch >= 0
    m = 2;
else
    m = 1;
end
%伪随机幅值
A = 1;
%PWM延迟拍数
Beat= 1;

opts = bodeoptions;
opts.FreqUnits = 'Hz';       
opts.PhaseWrapping = 'on';
%阶跃持续时间
Step = 5000;


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
Ud = [t_total, input_total];
Uq = [t_total, zeros(length(input_total),1)];

%仿真辨识
out = sim("Model.slx");
output = out.Id.Data(Step + 1:end) - mean(out.Id.Data(Step + 1:end));
output = squeeze(output);
[sys_tf,Gest] = DataProcess(input_prbs - DC,output,k,m,Ts);
Ld_o = -Rs*Ts./log((pole(sys_tf)));
fprintf('辨识结果：%.5f mH\n相对偏差：%.5f%%\n\n', ...
    Ld_o(1)*1000, (Ld_o(1)-Ld)/Ld*100);

%% 多工作点辨识
% Ld_data = zeros(length(-5:0.1:5),1);
% Error = zeros(length(-5:0.1:5),1);
% AF = zeros(length(-5:0.1:5),1);
% PF = zeros(length(-5:0.1:5),1);
% BFT = zeros(length(-5:0.1:5),1);
% point=1;
% U=-5:0.1:5;
% for i=-5:0.1:5
%     DC = i;
%     fprintf('工作点%.1fV：\n',DC);
%     %生成输入信号
%     input_step = ones(Step,1)*DC;
%     t_step = (0:Ts:(Step-1)*Ts)';
%     input_prbs = idinput(2^N-1,'prbs');
%     for j=1:length(input_prbs)
%         if(input_prbs(j)==1)
%             input_prbs(j)=A+DC;
%         else
%             input_prbs(j)=-A+DC;
%         end
%     end
%     t_prbs = (Step*Ts : Ts : (Step + 2^N - 2)*Ts)';
%     t_total = [t_step; t_prbs];
%     input_total = [input_step; input_prbs];
%     Ud = [t_total, input_total];
%     Uq = [t_total, zeros(length(input_total),1)];
%     %仿真辨识
%     out = sim("Ld_Identification_sim_simscape.slx");
%     output = out.Id.Data(Step + 1:end) - mean(out.Id.Data(Step + 1:end));
%     output = squeeze(output);
%     [sys_tf,~,r_AF,r_PF,r_BFT] = DataProcess(input_prbs - DC,output,k,m,Ts);
%     close all;
%     Ld_o = -Rs*Ts./log((pole(sys_tf)));
%     fprintf('辨识结果：%.5f mH\n相对偏差：%.5f%%\n\n', ...
%         Ld_o(1)*1000, (Ld_o(1)-Ld)/Ld*100);
%     Ld_data(point,1)=Ld_o(1);
%     Error(point,1)=(Ld_o(1)-Ld)/Ld*100;
%     AF(point,1)=r_AF*100;
%     PF(point,1)=r_PF*100;
%     BFT(point,1)=r_BFT*100;
%     point = point+1;
% end
% figure(Name='拟合程度');
% plot(U,AF,U,PF,U,BFT);
% title('拟合程度分析');
% xlabel('D轴工作点电压/V');
% ylabel('拟合程度/%');
% legend('幅频相关系数','相频相关系数','拟合优度');
% figure(Name='相对误差');
% plot(U,Error);
% title('不同工作点电感辨识相对误差');
% xlabel('D轴工作点电压/V');
% ylabel('相对误差/%');