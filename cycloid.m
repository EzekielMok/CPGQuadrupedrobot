%Ezekiel
%2023年3月25日
%% 清理变量
clc;
clear;
close all;
n=2000;
t_total = 6*pi/2;
t=linspace(0,t_total,n);%转一圈
ts=t_total/n;
 trace_sq = nan+t;%轨迹配置初始化，空值
  trace_xyz = nan+[t;t;t];%轨迹配置初始化，空值
    trace_abg= nan+[t;t;t];%轨迹配置初始化，空值
  trace_triang = nan+t;%轨迹配置初始化，空值
 k=0;
 omega = 2*pi;% 控制步态周期
 siggnal_sq_intergration1 = 0;
  siggnal_sq_intergration2 = 0;
 S = 50;
 H = 40;
 Tm = 0.5;
for j=t
    k=k+1;
    singnal_sq1_1 =0.5*square(omega*j)+0.5;
    singnal_sq1_2 =0.5*square(omega*j/2)+0.5;
    singnal_sq2_1 =0.5*square(2*omega*j)+0.5;
    singnal_sq2_2 =0.5*square(omega*j)+0.5;
    siggnal_sq_intergration1 = siggnal_sq_intergration1 + ts*logic_clock(singnal_sq1_1);
    singnal_triang1 = siggnal_sq_intergration1*logic_clock(singnal_sq1_2);
    siggnal_sq_intergration2 = siggnal_sq_intergration2 + ts*logic_clock(singnal_sq2_1);
    singnal_triang2 = siggnal_sq_intergration2*logic_clock(singnal_sq2_2);
    trace_triang(k)=  siggnal_sq_intergration1;
    x_3d = cycloid_x(singnal_triang1,Tm,S);
    z_3d = 150-cycloid_y(singnal_triang2,Tm,H);
   y_3d = 49;
    [alfa_dg, beta_dg, gamma_dg]=xyztoang( x_3d, y_3d,  z_3d);
    trace_xyz(:, k)= [x_3d;y_3d;z_3d];
    trace_abg(:, k)= [alfa_dg; beta_dg; gamma_dg];
    trace_sq(k)= singnal_sq1_1;
%  [alfa, beta, gamma]=xyztoang(x, y, z);
end
subplot(4,1,1)
plot(t,trace_sq(:),'red')
% title('Kimura振荡器输出')
ylabel('pulse')
axis([0,t_total ,-0.5,1.5])%XY坐标均衡
grid on;
subplot(4,1,2)
plot(t,trace_triang(:),'blue')
% title('Kimura振荡器输出')
ylabel('trangpulse')
axis([0,t_total ,0,1])%XY坐标均衡
grid on;
subplot(4,1,3)
% title('Kimura振荡器输出')
plot(trace_xyz(1, :),trace_xyz(3, :),'red')
ylabel('XY-plot')
axis([0,50 ,100,160])%XY坐标均衡
grid on;
subplot(4,1,4)
% title('Kimura振荡器输出')
plot(t, trace_abg(1, :), 'blue', t, trace_abg(2, :), 'red',t, trace_abg(3, :), 'green')
ylabel('theta')
axis([0,t_total ,-3,2])%XY坐标均衡
legend('alfa','beta','gamma'); 
grid on;
function y = logic_clock(u)
if u==0
    y = -1;
else
    y=1;
end
end
function x = cycloid_x(u,Tm,S)
    if u>0
        x =S*(u/Tm-1/(2*pi)*sin(2*pi*u/Tm));
    else 
        u=-u;
         x =S*(u/Tm-1/(2*pi)*sin(2*pi*u/Tm));
    end
end
function y = cycloid_y(u,Tm,H)
    n=4;
    if u>0
        y=2*H*(u/Tm-1/(n*pi)*sin(n*pi*u/Tm)).*(u>=0&u<Tm/2)+...
            2*H*(1-u/Tm+1/(n*pi)*sin(n*pi*u/Tm)).*(u>=Tm/2&u<Tm);
    else 
        y=0;
    end
end
function [alfa, beta, gamma]=xyztoang(x, y, z)
    h_up=49.0;
    h_mid=125.5899;
    h_low=116.0+20;
    dyz=sqrt(y.^2+z.^2);
    lyz=sqrt(dyz.^2-h_up.^2);
    gamma_yz=-atan(y./z);
    gamma_h_offset=-atan(h_up./lyz);
    gamma=gamma_yz-gamma_h_offset;
    lxzp=sqrt(lyz.^2+x.^2);
    n=(lxzp.^2-h_low^2-h_mid.^2)/(2*h_mid);
    beta=-acos(n/h_low);
    alfa_xzp=-atan(x/lyz);
    alfa_off=acos((h_mid+n)/lxzp);
    alfa=(alfa_xzp+alfa_off);
end
