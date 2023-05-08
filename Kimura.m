%Ezekiel
%2023年3月25日
%% 清理变量
% clc;
% clear;
% close all;
%% 参数设置
leg_num=4;%足的个数
mm=8;%m项反馈信息
Tr=0.04;%上升时间参数
Ta=0.38;
c=0.23;%直流输入
b=-2.0;%适应系数
a=-1.0;%抑制系数
phai=1;
%%
gait=2;
w_w=zeros(leg_num,leg_num);%定义步态矩阵
switch(gait)
    case 1
       w_w=[0 -1 1 -1;-1 0 -1 1;1 -1 0 -1;-1 1 -1 0];
    case 2
        w_w=[0 -1 -1 1;-1 0 1 -1;-1 1 0 -1;1 -1 -1 0];
    case 3
        w_w=[0 1 -1 -1;1 0 -1 -1;-1 1 0 -1;-1 -1 1 0];
    case 4
        w_w=[0 -1 -1 -1;-1 0 -1 -1;-1 -1 0 -1;-1 -1 -1 0];
    otherwise
        w_w=[0 -1 1 -1;-1 0 -1 1;1 -1 0 -1;-1 1 -1 0];   
end
%%
vf=zeros(leg_num,1);
ve=zeros(leg_num,1);
Vfe=[vf,ve];
%%
yf=zeros(leg_num,1);
ye=zeros(leg_num,1);
Yfe=[yf,ye];
Yef=[ye,yf];
Y=zeros(leg_num,1);%输出的关节映射矩阵
%%
uf=rand(leg_num,1)*c/10;
ue=rand(leg_num,1)*c/10;
Ufe=[uf,ue];
%%
hf=zeros(mm,1);
he=zeros(mm,1);
Hfe=[hf,he];
%%
ss=zeros(leg_num,mm);
WYf=zeros(leg_num,1);
WYe=zeros(leg_num,1);
SHf=zeros(leg_num,1);
SHe=zeros(leg_num,1);
%%
E=zeros(leg_num,2);
for i=1:1:leg_num
    E(i,1)=1;
    E(i,2)=1;
end
%% 关节参数，此处参数需实测，试验数据参考书
h=0.02;  %抬腿高度
v=1;   %行走速度
T=0.4; %步态周期
S=v*T; %步长
l=0.4;   %腿节长度
theta0=30/180*pi;%髋关节和膝关节平衡位置与垂直线夹角 
L=2*l*cos(theta0);%髋关节与足端之间长度
Ah=8.3;%asin((Beta*S)/(2*L));%髋关节摆动幅度
Ak=5.3;%acos((l*cos(theta0)-h)/l)-theta0;%膝关节摆动幅度

%% 时间，分辨率设置
%对时间进行算法积分
u0=[0.009 0.0021 0.0134 0.0022;
       0.0272 0.0003 0.0206 0.0057;
       0.0199 0.0067 0.0171 0.0111;
       0.0022 0.0151 0.0149 0.0081;];%初值会影响相序
n=2000;
Foot_end_x=zeros(n);%足端
Foot_end_y=zeros(n);
t_total = 8*pi/2;
t=linspace(0,t_total,n);%转一圈
ts=t_total/n;
% 积分计算
 uf = u0(:,1);
 vf = u0(:,2);
 ue = u0(:,3);
 ve = u0(:,4);
 shsum = zeros(leg_num,1);
 k = 0;
 trace_Y_Y = nan+[t;t;t;t];%轨迹配置初始化，空值
for j=t
    k=k+1;
    [WYf, WYe] = wye_sumeq(w_w, leg_num,  yf, ye);
    [uf_dot, vf_dot,ue_dot, ve_dot] = kimura_stateq(Tr, Ta, uf, vf, ue, ve, yf, ye,a, b, c,WYf, WYe, shsum);
    uf =uf + ts.*uf_dot;
    vf =vf + ts.*vf_dot;
    ue =ue + ts.*ue_dot;
    ve =ve + ts.*ve_dot;
    Y_Y = uf - ue;
    trace_Y_Y(:,k)= Y_Y ;
    [yf, ye] = returnyfe(uf, ue, leg_num);
end
yPlot = trace_Y_Y';
Judge1=diff(yPlot(:,1));
Judge2=diff(yPlot(:,2));
Judge3=diff(yPlot(:,3));
Judge4=diff(yPlot(:,4));
kplot=zeros(length(yPlot(:,1)),4);
thetah(1)=max(yPlot(:,1));
thetah(2)=max(yPlot(:,2));
thetah(3)=max(yPlot(:,3));
thetah(4)=max(yPlot(:,4));
u=2;%影响曲线幅值，幅值为开根号
for j=1:1:length(yPlot(:,1))-1
    if Judge1(j)>0
       kplot(j,1)=1*(Ak/Ah)*(1-(yPlot(j,1)/thetah(1))^2);
    else
       kplot(j,1)=0;
    end
    if Judge2(j)>0
       kplot(j,2)=1*(Ak/Ah)*(1-(yPlot(j,2)/thetah(2))^2);
    else
       kplot(j,2)=0;
    end
    if Judge3(j)>0
       kplot(j,3)=-1*(Ak/Ah)*(1-(yPlot(j,3)/thetah(3))^2);
    else
       kplot(j,3)=0;
    end
    if Judge4(j)>0
       kplot(j,4)=-1*(Ak/Ah)*(1-(yPlot(j,4)/thetah(4))^2);
    else
       kplot(j,4)=0;
    end
        Foot_end_x(j)=sin(kplot(j,4)/sqrt(u)*Ak*(pi/180)+theta0)*l-sin(theta0+yPlot(j,4)/sqrt(u)*Ah*(pi/180))*l;
        Foot_end_y(j)=L-(cos(yPlot(j,4)/sqrt(u)*Ah*(pi/180)+theta0)*l+cos(kplot(j,4)/sqrt(u)*Ak*(pi/180)+theta0)*l);
end
figure(1)
subplot(4,1,1)
plot(t,yPlot(:,1),t,kplot(:,1),'green')
% title('Kimura振荡器输出')
ylabel('LF')
axis([0,10,-2,2])%XY坐标均衡
legend('thetah','thetak');
grid on;
subplot(4,1,2)
plot(t,yPlot(:,2),t,kplot(:,2),'green')
ylabel('RF')
axis([0,10,-2,2])%XY坐标均衡
legend('thetah','thetak');
grid on;
subplot(4,1,3)
plot(t,yPlot(:,3),t,kplot(:,3),'green')
ylabel('RB')
axis([0,10,-2,2])%XY坐标均衡
 legend('thetah','thetak');
grid on;
subplot(4,1,4)
plot(t,yPlot(:,4),t,kplot(:,4),'green')
xlabel('时间（t/s）')
ylabel('LB')
axis([0,10,-2,2])%XY坐标均衡
legend('thetah','thetak');
grid on;

figure(2)
plot(Foot_end_x(1:120), Foot_end_y(1:120))
grid on
hold off
function [yf, ye] = returnyfe(uf, ue,legnum)
    yf = zeros(legnum,1);
    ye= zeros(legnum,1);
    for i= 1:1:legnum
        yf(i) = max(uf(i), 0);
        ye(i) = max(ue(i), 0);
    end
end
function [uf_dot, vf_dot,ue_dot, ve_dot] = kimura_stateq(Tr, Ta, uf, vf, ue, ve, yf, ye,a, b, c,wyf_sum, wye_sum, shsum)
    uf_dot =(1/Tr)*(-uf+b*vf+a*ye+wyf_sum+shsum+c);
    vf_dot = (1/Ta)*(yf-vf);
    ue_dot = (1/Tr)*(-ue+b*ve+a*yf+wye_sum+shsum+c);
    ve_dot =  (1/Ta)*(ye-ve);
end
function [WYf, WYe] = wye_sumeq(w, legnum,  yfaugment, yeaugment)
WYf= zeros(legnum,1);
WYe= zeros(legnum,1);
for i_i =1:1:legnum
    for j_j=1:1:legnum
        WYf(i_i)=WYf(i_i)+w(i_i, j_j)*yfaugment( j_j);
        WYe(i_i)=WYe(i_i)+w(i_i, j_j)*yeaugment( j_j);
    end
end
end
function shsum = shsumeq(s, legnum, haugment)
shsum = zeros(legnum,1);
for i =1:1:legnum
    for j=1:1:legnum
        shsum(i)=shsum(i)+s(i,j)*haugment(j);
    end
end
end 