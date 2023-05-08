clc;clear;close all;
lx = 1; lr = 0.1;  % 连杆的长度和半径
gy = 9.81;  % 重力加速度(这里在y轴方向，因为z轴是关节转动轴)
fvis = 0; fcou = 0; % 粘性摩擦系数和库伦摩擦系数
rod = Cuboid([lx,lr,lr]); % 建立成长方体
Irod = rod.inertia; % 可以直接获得转动惯量
rod = Cylinder(lr,lx);
Rcyl = SO3.rpy([0 -pi/2 0]);
Irod = Rcyl.R*rod.inertia*Rcyl.R';
% 设定dh参数，质量，关节质心距离，关节约束，转动惯量，摩擦系数
dpm = {'a', lx, 'm', rod.mass, 'r', [-lx/2,0,0],...
 'qlim', [-pi/2, pi/2],'I', Irod,...
    'B', fvis, 'Tc', [fcou -fcou]};
% 建立二连杆机械臂，设定回转关节和重力
r = SerialLink([Revolute(dpm{:}),Revolute(dpm{:})],...
'name','two-link','gravity',[0 gy 0]);
ws = [-4 4 -4 4 -4 4]; % 设置工作空间
plotopt = {'workspace', ws, 'nobase', 'notiles', 'noshading', 'noshadow', 'nowrist','top'}; % 设置绘图参数

% 给定初始关节角位置和速度
qz = zeros(1,2);
qd = zeros(1,2);
h = r.plot(qz,plotopt{:});
% 设定仿真时间10s
y0 = [qz,qd]'; tspan = [0 10]; 
% 刚体仿真用ode15s比较快，用accel函数直接获得关节加速度
tic
[tlist,ylist] = ode15s(@(t,y) [y(r.n+1:end);r.accel(y(1:r.n)',y(r.n+1:end)',zeros(1,r.n))],tspan,y0);
toc
% 差不多6s可出结果，然后画出机械臂即可
ws = [-4 4 -4 4 -4 4];
plotopt = {'workspace', ws, 'nobase', 'notiles', 'noshading', 'noshadow', 'nowrist','top'};
h = r.plot(ylist(:,1:r.n),plotopt{:});
