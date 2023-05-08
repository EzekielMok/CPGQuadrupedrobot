%Ezekiel
%2023年3月25日
%% 清理变量
clc;
clear;
close all;
n=2000;
r =1;
for index = 0:1:n
    theta = 5*pi/1000*index;
    x(index+1) = r*(theta - sin(theta));
    y(index+1) = r*(1 - cos(theta));
end
plot (x, y,'-r','linewidth',1);
axis equal
grid on
xlabel('X')
ylabel('Y')