% 时域信号分析函数
function [r,BFT] = Signal_Analyse (Signal_Rec, Signal_Sim)
%输入参数定义
%Signal_Rec    实测信号
%Signal_Sim    仿真信号
%输出参数定义
%r             相关系数
%BFT           拟合优度

cov_sim_rec = cov(Signal_Rec, Signal_Sim);
cov_sim_rec=cov_sim_rec(1,2);
r=cov_sim_rec/(sqrt(var(Signal_Rec))* ...
    sqrt(var(Signal_Sim)));
fprintf('相关系数：%.5f%%\n',r*100);
BFT = (1 - norm(Signal_Rec - Signal_Sim) / ...
    norm(Signal_Rec - mean(Signal_Rec)));
fprintf('拟合优度：%.5f%%\n',BFT*100);

end
