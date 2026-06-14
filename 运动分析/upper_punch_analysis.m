%% 15吨压片机 上冲头机构 运动分析（解析法）
% 机构链: 凸轮 -> 下摆杆(支点O2) -> 长连杆L1 -> 上摆杆(支点O3,直杆杠杆)
%         -> 小连杆L2 -> 冲头(销E, 沿 X=0 竖直导向)
% 输入: 主轴转角 phi(0~360°), 凸轮使下摆杆作摆线摆动 psi(phi)
% 输出: 冲头尖端 位移 s、速度 v、加速度 a 随 phi(=时间) 的变化曲线
% 数据取自重建方案 结果.json (单位 mm / 度)；与该数据已数值核对(行程=95.000)
clear; clc; close all;

%% 1. 机构参数
O2   = [100, 35];        % 下摆杆支点
LB   = 150;             % 下摆杆 连杆臂 (O2->B)
lam0 = -9.0236;         % 装填位(phi=0) 连杆臂方位角 (deg, 自+X逆时针)
psimax = 18.2176;       % 下摆杆 最大摆角 (deg)
O3   = [140, 457.039];  % 上摆杆支点
R3   = 70;              % 上摆杆 短臂 (O3->D)
R3p  = 140;             % 上摆杆 长臂 (O3->C, 与短臂反向共线)
L1   = 423.930;         % 长连杆 (B->D)
L2   = 75;              % 小连杆 (C->E)
Lp   = 130;             % 冲头销E -> 尖端 长度
n    = 25;              % 主轴转速 r/min
w    = n*2*pi/60;       % 主轴角速度 rad/s (=2.618)

%% 2. 凸轮运动规律 f(phi) —— 摆线
% 相位: 0~120 近休(高位) | 120~200 推程(下压) | 200~240 远休(保压)
%       240~300 回程(退回) | 300~360 近休
phi = linspace(0,360,3601);          % deg
cyc = @(tau) tau - sin(2*pi*tau)/(2*pi);
f = zeros(size(phi));
for k = 1:numel(phi)
    p = phi(k);
    if     p>=120 && p<200, f(k) = cyc((p-120)/80);
    elseif p>=200 && p<240, f(k) = 1;
    elseif p>=240 && p<300, f(k) = 1 - cyc((p-240)/60);
    else,                    f(k) = 0;
    end
end
psi = psimax * f;                    % 下摆杆摆角 (deg)

%% 3. 位置分析（逐点解机构闭环）
tip = zeros(size(phi));  th3 = zeros(size(phi));  Ddet = [];
for k = 1:numel(phi)
    lam = deg2rad(lam0 + psi(k));
    B   = O2 + LB*[cos(lam), sin(lam)];          % 下摆杆连杆臂端
    % D: 圆(O3,R3) 与 圆(B,L1) 的交点
    dvec = B - O3;  dd = norm(dvec);
    a_ = (R3^2 - L1^2 + dd^2)/(2*dd);
    h  = sqrt(max(R3^2 - a_^2, 0));
    M  = O3 + a_*dvec/dd;  perp = [-dvec(2), dvec(1)]/dd;
    D1 = M + h*perp;  D2 = M - h*perp;
    if isempty(Ddet)                              % 用装填位已知解选根
        ref = [205.85, 433.29];
    else
        ref = Ddet;                               % 连续性选根
    end
    if norm(D1-ref) < norm(D2-ref), D = D1; else, D = D2; end
    Ddet = D;
    th = atan2(D(2)-O3(2), D(1)-O3(1));  th3(k) = rad2deg(th);
    C  = O3 - R3p*[cos(th), sin(th)];             % 上摆杆长臂端(反向)
    Ey = C(2) - sqrt(L2^2 - C(1)^2);              % 冲头销E 在 X=0, 取下解
    tip(k) = Ey - Lp;                             % 冲头尖端 Y
end

%% 4. 速度、加速度（对时间求导, phi = w*t 均匀）
t = deg2rad(phi)/w;                  % 时间 s
v = gradient(tip, t);                % 速度 mm/s
a = gradient(v,   t);                % 加速度 mm/s^2

%% 5. 关键结果
stroke = max(tip) - min(tip);
[vmax,iv] = max(abs(v));  [amax,ia] = max(abs(a));
fprintf('冲头行程 = %.3f mm  (设计值 95)\n', stroke);
fprintf('最大速度 |v|max = %.1f mm/s  @ phi = %.0f°\n', vmax, phi(iv));
fprintf('最大加速度 |a|max = %.0f mm/s^2  @ phi = %.0f°\n', amax, phi(ia));

%% 6. 绘制 位移-速度-加速度 曲线
figure('Color','w','Position',[100 80 720 820]);
subplot(3,1,1); plot(phi,tip,'b-','LineWidth',1.6); grid on; box on;
ylabel('位移 s /mm'); title('上冲头尖端 运动分析（解析法）');
xlim([0 360]); xticks(0:60:360);
subplot(3,1,2); plot(phi,v,'r-','LineWidth',1.6); grid on; box on;
ylabel('速度 v /(mm·s^{-1})'); xlim([0 360]); xticks(0:60:360); yline(0,'k:');
subplot(3,1,3); plot(phi,a,'m-','LineWidth',1.6); grid on; box on;
ylabel('加速度 a /(mm·s^{-2})'); xlabel('主轴转角 \phi /(°)');
xlim([0 360]); xticks(0:60:360); yline(0,'k:');

%% 7. 导出数据（供说明书/与SolidWorks Motion对比）
T = table(phi.', psi.', tip.', v.', a.', ...
    'VariableNames', {'phi_deg','psi_deg','s_mm','v_mmps','a_mmps2'});
writetable(T, '上冲头运动分析数据.csv');
fprintf('已导出 上冲头运动分析数据.csv\n');
